classdef splitCoder < dynamicprops
  %
  % Container class for generating c-code solvers using matlab coder
  %
  
  properties
    prob   % splitProb to generate
    
    % Parameters common to all methods
    A,b,pB,Q,f,pF,c,L,l,pL,pQ,pf
    numDefaultParams = 12
    proxInd
    
    % Optimization parameters common to all methods
    maxItr    = 1e3
    dualTol   = 1e-3
    primalTol = 1e-3
    
    % Name of file to generate
    fname
    pathstr
    
    % Template name
    templateName
    
    % Output arguments
    outArgs
    
    % Parametric input arguments
    inArgs
    
    % Parameters that are available to the generation function, and to the compiled code
    inConstArgs
    
    % Parameters that are available to the generation function, and to the compiled code
    compileTimeParameteres
    
    % File handle to the generated file
    fHdl
  end

  % Private properties
  properties
    functions % Helper functions to be added at the end of the generated code
              % Elements are all strings:
              % functions.header 
              % functions.body
              % functions.comment 
              
    allocVars % Internal variables that require memory allocation
              %   .name
              %   .m 
              %   .n
              
    headerFuncs % List of functions to execute in order to 
                % generate code at the point >>> HEADER >>> after a
                % first-run of the compiler
  end

  
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function cdr = splitCoder(prob, templateName, fname)
      % Set the default arguments
      argNames = {'f','b','l','A','L','Q','pB','pF','pL','dualTol','primalTol','maxItr'};
      argDescs = {'','','','','','','','','','Dual tolerance','Primal tolerance','Maximum number of iterations'};
      for i = 1:length(argNames)
        cdr.inConstArgs(i).name = argNames{i};
        cdr.inConstArgs(i).desc = argDescs{i};
      end
      
      argNames = {'x_'};
      argDescs = {'Optimizer'};
      for i = 1:length(argNames)
        cdr.outArgs(i).name = argNames{i};
        cdr.outArgs(i).desc = argDescs{i};
      end
      
      % Process the problem
      cdr.prob  = prob;
      cdr.genData;
      
      [cdr.pathstr,cdr.fname,~] = fileparts(fname);
      if isempty(cdr.pathstr)
        cdr.pathstr = '.';
      end
      
      % Convert the name into a path + filename
      [pathstr,~,~]=fileparts(mfilename('fullpath'));
      cdr.templateName = sprintf('%s%stemplates%stemplate_%s.m', pathstr, filesep, filesep, templateName)
      
      % Open the template file and parse the IO block
      template = cdr.loadTemplate;
      
      % Find row where IO block ends
      indIO     = find(~cellfun('isempty',regexpi(template,'\s*%+\s*>+\s*IO\s*>+.*')));
      assert(length(indIO) == 1, 'Invalid template file: Could not find IO block')
      
      % Run the IO block
      cdr.executeBlock({template{1:indIO-1}});
    end
  end

  methods (Access = private)
    addVar(cdr, name, m, n) % Define a new variable
    allocateMemory(cdr)     % Allocate memory for internal variables
  end
  
  methods (Access = private)
    function genFunctionHeader(cdr)
      % Generate the function header
      cdr.print('function [%s] = %s(%s)%%#codegen\n', cdr.csvList({cdr.outArgs.name}), ...
        cdr.fname, ...
        cdr.csvList({cdr.inArgs.name,cdr.inConstArgs.name}));

      % Add functions that will be called after generation and whose output
      % will be placed at the >>> HEADER >>> block
      
      % Add description of the file being generated
      cdr.addHeaderFunc(@() cdr.functionDescription(cdr.fHdl, '%'));
      cdr.addHeaderFunc(@() cdr.print('\n'));
      
      % Add code for timing
      cdr.addHeaderFunc(@()...
        cdr.print(['if coder.target(''MEX''), coder.ceval(''split_tic''); else tic; end\n']));
      cdr.addHeaderFunc(@() cdr.print('\n'));
      
      % Add code to allocate memory 
      cdr.addHeaderFunc(@()cdr.allocateMemory);
      cdr.addHeaderFunc(@() cdr.print('\n'));
    end
    
    function genFunctionEnd(cdr)
      cdr.print('tm = 0.0;\n');
      cdr.print('if coder.target(''MEX'')\n')
      cdr.print('  tm = coder.ceval(''split_toc'');\n');
      cdr.print('else\n')
      cdr.print('  tm = toc;\n')
      cdr.print('end\n')
      cdr.print('end\n')
      
      for i = 1:length(cdr.functions)
        cdr.print('\n\n')
        cdr.print(cdr.functions(i).header)
        cdr.print('\n')
        for j = 1:length(cdr.functions(i).body)
          cdr.print(cdr.functions(i).body{j})
          cdr.print('\n')
        end
        cdr.print('end\n')
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function executeBlock(cdr, block)
      % Create a matlab file, copy the code in and execute it
      name = 'tmp.m';
      f = fopen(name,'w+');
      for i = 1:length(block)
        fprintf(f,strrep([block{i} '\n'],'%','%%'));
      end
      fclose(f);
      
      % Evaluate the code as a script
      rehash % Need to re-search the path, or matlab won't know about the new tmp file here
      tmp
      delete(name);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function addConst(cdr, input, desc)
      % Add a constant that will be available in the generated code
      I = find(strcmp({cdr.inConstArgs.name},input));
      if isempty(I)
        cdr.addprop(input);
        I = length(cdr.inConstArgs)+1;
      end        
      cdr.inConstArgs(I).name = input;
      cdr.inConstArgs(I).desc = desc;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function addParameter(cdr, input, desc)
      % Add a constant that must be specifid at compile-time
      I = length(cdr.compileTimeParameteres)+1;
      if ~isempty(cdr.compileTimeParameteres)
        I = find(strcmp({cdr.compileTimeParameteres.name},input));
        if isempty(I)
          cdr.addprop(input);
          I = length(cdr.compileTimeParameteres)+1;
        end
      else
        cdr.addprop(input);
        I = 1;
      end
      cdr.compileTimeParameteres(I).name = input;
      cdr.compileTimeParameteres(I).desc = desc;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function addInput(cdr, name, sz, desc)
      % Add an input to the generated function
      I = find(strcmp({cdr.inConstArgs.name},name));
      if isempty(I)
        I = length(cdr.inArgs)+1;
      end
      cdr.inArgs(I).name = name;
      cdr.inArgs(I).size = sz;
      cdr.inArgs(I).desc = desc;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function addOutput(cdr, name, desc)
      % Final output in the generated file
      I = find(strcmp({cdr.outArgs.name},name));
      if isempty(I), I = length(cdr.outArgs)+1; end
      cdr.outArgs(I).name = name;
      cdr.outArgs(I).desc = desc;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function addFunction(cdr, header, body, comment)
      % Add a function that will be added to the end of the file
      if nargin < 4, comment = ''; end
      ind = length(cdr.functions) + 1;
      cdr.functions(ind).header  = header;
      cdr.functions(ind).body    = body;
      cdr.functions(ind).comment = comment;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function addHeaderFunc(cdr, func)
      % Add a function that will be evaluated after the code has been
      % generated. Any output from this function will be placed at the mark
      % >>>> HEADER >>>>
      % func must be a function handle that takes no arguments
      cdr.headerFuncs{end+1} = func;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function mname = getFileName(cdr)
      pth = cdr.pathstr;
      if ~exist(pth, 'dir')
        mkdir(pth);
      end
        
      mname = [pth filesep cdr.fname '.m'];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function openGenFile(cdr, doClear, mode)
      % Create the m-file for output
      mname = cdr.getFileName();
      perm = 'a+';
      if doClear, perm = 'w+'; end
      if nargin == 3
        f = fopen(mname, mode)
      else
        f = fopen(mname, perm);
      end
      assert(f > 0, 'Could not open file %s', cdr.fname)
      cdr.fHdl = f;
    end
    function closeGenFile(cdr)
      fclose(cdr.fHdl);
      cdr.fHdl = [];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function print(cdr, varargin)
      if length(varargin) == 1 && ~isempty(regexp(varargin{1},'\s*%+.*'))
        % This is a comment - escape the %'s
        fprintf(cdr.fHdl, strrep(varargin{1},'%','%%'));
      else
        fprintf(cdr.fHdl, varargin{:});
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function c = loadTemplate(cdr)
      t = fopen(cdr.templateName, 'r'); assert(t > 0, 'Could not open template file')
      c = textscan(t, '%s', 'delimiter','\n','whitespace','');
      c = c{1};
      fclose(t);
    end
  end
  
  methods (Access = private)
    function cdr = genData(cdr)
      %----------------------------------------------------------------
      % Generate the data required for code-export
      %----------------------------------------------------------------
      
      f = fields(cdr.prob.dat);
      for i = 1:length(f)
        eval(sprintf('cdr.%s = full(cdr.prob.dat.%s);\n', f{i}, f{i}));
      end
      cdr.f = cdr.f(:);
      cdr.pF = cdr.pF';
      
      dat = cdr.prob.dat;
      if ~isfield(dat, 'pB')
        cdr.pB = zeros(size(dat.A,1),0);
        cdr.pF = zeros(size(dat.A,2),0);
      end
      
      prox    = cdr.prob.prox;
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
      end
      L = full(L); l = full(l);
      if ~isfield(dat, 'pB')
        pL = zeros(size(L,1),0);
      end
      
      cdr.L  = L;
      cdr.l  = l;
      cdr.pL = pL;
      cdr.proxInd = proxInd;
      
      if size(pL,2) > 0
        cdr.addInput('par', [size(pL,2),1], 'Parameter variable');
      end
    end
    
    function str = genProxFunction(cdr, x, y, rho, conj)
      %----------------------------------------------------------------
      % genProxFunction(cdr, x, y, rho, conj);
      %
      % Generate code to evaluate the proximal functions
      %
      % x, y and rho are strings containing the variable names
      %    y: output
      %    x: input
      %  rho: weight
      % conj: [false] if true, generate the proximal functions for the
      % conjugate functions
      %----------------------------------------------------------------
      if nargin < 5, conj = false; end
      
      prox = cdr.prob.prox;
      j = 1;
      for i = 1:length(prox)
        I = sprintf('%i:%i', cdr.proxInd(i)+1,cdr.proxInd(i+1));
        if conj, func = @(a,b)prox(i).conjStr(a,b);
        else     func = @(a,b)prox(i).funcStr(a,b);
        end
        str{j} = sprintf('%s(%s) = %s;', y, I, func([x '(' I ')'],rho));
        j = j + 1;
      end
    end
  end
  
  methods (Static)
    function str = csvList(varargin)
      %----------------------------------------------------------------
      % Return a string of comma seperated arguments
      %----------------------------------------------------------------
      list = {};
      for i = 1:length(varargin)
        list = {list{:} varargin{i}{:}};
      end
      % Print comma-seperated list
      str = sprintf('%s, ', list{1:end-1});
      str = [str list{end}];
    end
  end
  
  methods (Access = private)
    function genSparseLU(cdr,L,U,p,q,x,y,b)
      %
      % genSparseLU(L,U,p,q,x,y,b)
      %
      % Generate a sparse LU factorization solve
      %
      %  L,U,p,q - output from LU
      %  x,y,b   - [str] names of the variables involved
      %    x - output
      %    y - tmp variable
      %    b - rhs
      %
      
      n = size(L,1);
      
      assert(all(abs(diag(U))>0) && all(abs(diag(L))>0), 'Matrix is not full rank')
      
      % L \ b
      % str{end+1} = sprintf('y = coder.nullcopy(zeros(%i,1));', n);
      for i = 1:n
        % y(i) = (b(p(i))-L(i,1:i-1)*y(1:i-1)) / L(i,i);
        [~,I,V] = find(L(i,1:i-1));
        cdr.print('%s(%i) = %s%s;\n', y,i, ...
          cdr.genSparseMatVec(1/L(i,i),b,'',p(i),[],true),...
          cdr.genSparseMatVec(V/L(i,i), y, '', I,[],true));
      end
      
      % x(q) = U \ y
      for i = n:-1:1
        % x(q(i),1) = (y(i) - U(i,i+1:end)*x(q(i+1:end))) / U(i,i);
        [~,I,V] = find(U(i,i+1:end));
        J = q(i+1:end);
        cdr.print('%s(%i) = %s%s;\n', x, q(i), ...
          cdr.genSparseMatVec(1/U(i,i),y,'',i,[],true),...
          cdr.genSparseMatVec(-V/U(i,i),x,'',J(I),[],true));
      end
    end
    
    function genComputeParametric(cdr, par)
      % genComputeParametric(cdr, par)
      %
      % Generate code to compute the static problem from the parametric
      % parameter
      %
      
      if nnz(cdr.pF) > 0, cdr.print('f = pF*%s + f;\n', par); end
      if nnz(cdr.pL) > 0, cdr.print('l = pL*%s + l;\n', par); end
      if nnz(cdr.pB) > 0, cdr.print('b = pB*%s + b;\n', par); end
    end
  end
end