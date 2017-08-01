classdef pyCoderData < handle
  % Contains the setup data for the python code generation module
  properties
    % Arguments of the generated function (must be vectors)
    args = struct
    
    % Data available to the generator at compile time
    data = struct
    
    % Constants that will be generated in the code
    constants = struct
    
    % Settings passed to the generator
    settings = struct
  end
  
  
  methods
    function cdr = pyCoderData(prob)
      cdr.parseProb(prob); % Parse basic data
      
      % Add default values
      cdr.addConstant('MAXITR', 1e3);
      cdr.addConstant('DUALTOL', 1e-3);
      cdr.addConstant('PRIMALTOL', 1e-3);
      cdr.addConstant('ITR_PER_CONV_TEST', 10);
      cdr.addConstant('rho', 1.0);
      
      cdr.addSetting('fname', 'mySolver');
      cdr.addSetting('functionName', 'solve');
      
      cdr.addArg('par', size(cdr.data.pB,2));
    end
    
    function addData(cdr, name, value),     cdr.add('data',name,value); end
    function addArg(cdr, name, value),      cdr.add('args',name,value); end
    function addConstant(cdr, name, value), cdr.add('constants',name,value); end
    function addSetting(cdr, name, value),  cdr.add('settings',name,value); end
    
    function rmData(cdr, name, value),     cdr.rm('data',name,value); end
    function rmArg(cdr, name, value),      cdr.rm('args',name,value); end
    function rmConstant(cdr, name, value), cdr.rm('constants',name,value); end
    function rmSetting(cdr, name, value),  cdr.rm('settings',name,value); end
    
    function genData(cdr, fname)
      % writegenData(fname)
      %
      % Write the data to the file fname
      %
      
      if nargin < 2, fname = 'genData.py'; end
      
      f = fopen(fname,'w+');
      assert(f ~= -1, 'Could not open file %s', fname)
      
      fprintf(f, 'import numpy as np\n\n\n');

      types = {'args','constants','settings','data'};
      names = {'args','constants','settings','dat'};
      for i = 1:length(types)
        fprintf(f, '%s = ', names{i});
        cdr.printStruct(f, cdr.(types{i}));
        fprintf(f, '\n\n');
      end
      
      fclose(f);
    end
  end
  
  methods(Hidden=true, Access=private)
    function add(cdr, type, name, value)
      for str = {'args','data','constants','settings'}
        cdr.rm(str{:}, name);
      end
      cdr.(type).(name) = value;
    end
    
    function rm(cdr, type, name)
      if isfield(cdr.(type), name)
        fprintf('Removing %s from %s\n', name, type);
        cdr.(type) = rmfield(cdr.(type), name);
      end
    end
    
    function parseProb(cdr, prob)
      dat = prob.dat;
      
      if ~any(isfield(dat, {'pF','pB','pL'}))
        error('Can only generate solvers for Argetric problems')
      end
      
      % Fix matrix sizes
      dat.f  = dat.f(:);
      dat.pF = dat.pF';
      
      % Add empty Argetric matrices if they're not present
      if ~isfield(dat, 'pB')
        dat.pB = zeros(size(dat.A,1),0);
        dat.pF = zeros(size(dat.A,2),0);
      end
      
      % Build the prox functions
      prox    = prob.prox;
      nProx = length(prox);
      L = []; l = []; pL = []; w = [];
      proxInd = 0;
      for i = 1:nProx
        L  = [L;prox(i).L];
        l  = [l;prox(i).l];
        proxInd(i+1) = size(L,1);
        if isfield(dat, 'pB')
          pL = [pL;full(prox(i).pL)];
        end
        dat.prox(i).type = prox(i).typeName;
        dat.prox(i).ind  = proxInd(i);
        dat.prox(i).len  = size(prox(i).L,1);
        dat.prox(i).dat  = prox(i).dat;
      end
      L = full(L); l = full(l);
      if ~isfield(dat, 'pB')
        pL = zeros(size(L,1),0);
      end
      
      dat.L  = L;
      dat.l  = l;
      dat.pL = pL;
      dat.proxInd = proxInd;
      
      cdr.data = dat;
    end
    
    function printStruct(cdr, f, dat)
      fprintf(f, '{\n');
      names = fieldnames(dat);
      for i = 1:length(names)
        fprintf(f, '\t"%s": ', names{i});
        d = getfield(dat, names{i});
        if iscell(d)
          fprintf(f, '(');
          for k = 1:length(d)
            assert(isstr(d{k}), 'Can only handle cell arrays of strings')
            fprintf(f, '"%s", ', d{k});
          end
          fprintf(f, ')');
        elseif isnumeric(d)
          if isscalar(d) && isinf(d)
            fprintf(f, '"Inf"');
          else
            cdr.printMatrix(f, d);
          end
        elseif isstr(d)
          fprintf(f, '"%s"', d);
        elseif isstruct(d)
          if length(d) > 1
            fprintf(f, '(\n');
            for k = 1:length(d)
              cdr.printStruct(f, d(k));
              fprintf(f, ',');
            end
            fprintf(f, ')');
          else
            cdr.printStruct(f, d);
          end
        else
          error('Unknown field type %s\n', names{i})
        end
        fprintf(f,',\n');
      end
      fprintf(f, '}\n');
    end
    
    function printMatrix(cdr, f, M)
      if prod(size(M)) == 1 % Scalar
        fprintf(f, '%g', full(M));
      elseif prod(size(M)) == 0 % Empty - replicate the size of the zero matrix
        fprintf(f, 'np.zeros((%i,%i))', size(M,1), size(M,2));
      else
        fprintf(f, 'np.array([');
        for i = 1:size(M,1)
          fprintf(f, '[');
          for j = 1:size(M,2)
            fprintf(f, '%g, ', full(M(i,j)));
          end
          fprintf(f, '], ');
        end
        fprintf(f, '])');
      end
    end
  end
end
