classdef splitData < handle
  %
  % Class to write double and integer vectors to a binary data file
  %
  properties
    data  = struct([])
    
    fHdr  = -1;
    fname  = '';
  end
  
  methods
    function obj = splitData(fname)
      % Format: splitData(fname)
      %
      % fname : name of data file [default = splitData]
      % Creates two files: 
      %   - fname.dat = binary file containing the data
      %   - fname.h   = header file containing static data allocations for all the variables   
      
      obj.fHdr = fopen(sprintf('%s.h', fname), 'w+');
      obj.fname = fname;
      
      pl(obj.fHdr, '#ifndef __%s_h__', fname)
      pl(obj.fHdr, '#define __%s_h__', fname)
      
      pl(obj.fHdr, '#include <stdio.h>')
      pl(obj.fHdr, '#include <string.h>')
    end

    function add(obj, name, x, type, genStatic)
      % Add variable to be written to the data file
      
      % Check if variable already exists
      for i = 1:length(obj.data)
        if strcmp(obj.data(i).name, name)
          % Confirm that the data is the same
          if norm(obj.data(i).x - x) > 0
            error('Attempt to add variable %s to the data file twice with different data', name)
          end
          return % This variable already exists
        end
      end
      
      i = length(obj.data)+1;
      obj.data(i).name = name;
      obj.data(i).x    = x(:);
      if nargin < 4 || strcmpi(type,'double')
        obj.data(i).type = 0;
      else
        obj.data(i).type = 1;
      end
      
      if nargin < 5 || genStatic
        if obj.data(i).type == 0
          pl(obj.fHdr, 'double %s[%i];', name, length(x));
        else
          pl(obj.fHdr, 'int %s[%i];', name, length(x));
        end
        pl(obj.fHdr, '#define %s_len %i', name, length(x));
      end
    end
    
    function define(obj, name, def, type)
      % Define an item in the header file
      %
      % define(obj, name, def, type)
      %  #define name def
      %  
      % If type is unspecified, it will try to guess
      
      if nargin < 4
        if ischar(def)
          type = 'char';
        elseif isnumeric(def)
          if ceil(def) == def
            type = 'int';
          else
            type = 'float';
          end
        end
      end
          
      switch type
        case 'char'
          pl(obj.fHdr, '#define %s "%s"', name, def);
        case 'int'
          pl(obj.fHdr, '#define %s %i', name, def);
        case 'float'
          pl(obj.fHdr, '#define %s %f', name, def);
        otherwise
          error('Unknown type')
      end
    end
    
    function hdr = writeFile(obj)
      f = fopen(sprintf('%s.dat', obj.fname), 'w+');
      
      % Write the header
      % struct Header {
      %   unsigned int numVectors;
      %   char info[255];          // Text string describing the file
      % };
      hdr = sprintf(['Optimization data created by SPLIT\n'...
        'For more information on the SPLIT toolbox visit la.epfl.ch']);
      hdr = {length(obj.data), splitData.format_c_string(hdr, 255)};
      struct_write(f, 'Ic', hdr);
      
      % Write each variable
      for i = 1:length(obj.data)
        % struct Vector {
        %   unsigned int type; // 1 => int, 0 => double
        %   unsigned int len;  // Number of elements
        %   char name[100];     // Name of the variable
        %   void *data;      // Vector of sizeof(type)*len
        % };
        % typedef struct Vector Vector;
        d = {obj.data(i).type, ...
          length(obj.data(i).x),...
          splitData.format_c_string(obj.data(i).name,100),...
          full(obj.data(i).x)};
        if obj.data(i).type == 0
          struct_write(f, 'IIcd', d);
        elseif obj.data(i).type == 1
          struct_write(f, 'IIci', d);
        else
          error('Unknown data type')
        end
      end
      
      fclose(f);
      
      
      pl(obj.fHdr);
      pl(obj.fHdr);
      obj.loadData; % Write a function to load the data
      pl(obj.fHdr);
      pl(obj.fHdr);
      pl(obj.fHdr, '#endif');
      fclose(obj.fHdr); % Close the header file too
    end
    
    function loadData(obj)
      % Write out a function load all data into the memory
      pl(obj.fHdr, 'void loadData() {');
      pl(obj.fHdr, ' Data *dat = splitLoad("%s.dat");', obj.fname);
      for i = 1:length(obj.data)
        pl(obj.fHdr, ' copyVar(%s, dat, "%s");', obj.data(i).name, obj.data(i).name);
      end
      pl(obj.fHdr, 'freeData(dat);')
      pl(obj.fHdr, '}')
    end
    
    
    function [len,I,vec] = writeSparseMatrix(dat, A, Astr)
      %
      % Write out a constant sparse matrix rowwise
      %
      
      [n,m] = size(A);
      
      vec = [];
      I   = [];
      len = [];
      for i = 1:n
        J = find(A(i,:));
        a = full(A(i,J));
        
        len = [len length(J)];
        I   = [I J];
        vec = [vec a];
      end
      
      dat.add(sprintf('%s_nz_per_row',Astr), len, 'unsigned int');
      dat.add(sprintf('%s_ind',Astr), I-1, 'unsigned int');
      dat.add(sprintf('%s_dat',Astr), vec);
    end
  end
  
  methods (Static)
    function str = format_c_string(str, len)
      i = length(str);
      str = sprintf('%-*s',len, str);
      if i < len
        str(i+1) = 0;
      else
        str(len+1) = 0;
      end
    end
  end
end
