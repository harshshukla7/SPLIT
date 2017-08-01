classdef coderFunc < handle
  % 
  % Stores the description of a generated c-functions
  %
  
  properties
    prototype = '' % Function prototype
    body      = '' % Body of the function
    desc      = '' % Description of the function
    
    data      = struct('name',{},'x',{},'args',{}); % Stored data for this function
    functions = {}; % List of internal functions that need to be generated 
  end
  
  methods
    function func = coderFunc(varargin)
      % func = coderFunc(varargin)
      %
      % Defines a new function
      % 
      % Argument: sprintf-like text that gives the function prototype
      %
      func.prototype = sprintf(varargin{:});
    end
    
    %% These functions are used to generate the header and c-file
    function str = print_hfile(func, end_with_semicolon)
      % Returns string to go in the header file for this function
      str = '';
      if ~isempty(func.desc)
        str = sprintf('/*\n%s\n*/\n', func.desc);
      end
      str = sprintf('%s%s', str, func.prototype);
      if nargin < 2 || end_with_semicolon == true
        str = sprintf('%s;\n', str);
      end
    end

    function str = print_cfile(func)
      % Returns string to go in the c-file defining this function
      str = sprintf('%s {\n%s}\n\n', func.print_hfile(false), func.body);
    end

    function p(func, varargin)
      % Print to the body of the function
      if isempty(varargin), return; end
      func.body = [func.body sprintf(varargin{:})];
    end
    function pl(func, varargin)
      % Print a line to the body of the function
      if ~isempty(varargin)
        func.p(varargin{:})
      end
      func.p('\n')
    end
    
    function funcs = getFunctions(f)
      % funcs = getFunctions(f)
      %
      % Recursively returns a flat list of all sub-functions, including
      % this one
      %
      funcs = {f};
      for i = 1:length(f.functions)
        subFuncs = f.functions{i}.getFunctions;
        funcs = {funcs{:} subFuncs{:}};
      end
    end
    
    function add_func(f, func)
      % Add a sub-function that this function depends on
      f.functions{end+1} = func;
    end
  end
  
  methods(Access = 'protected', Hidden = true)
    function add_var(f, name, x, varargin)
      % Add a variable that needs to be stored 
      % This function just records the arguments, and passes them onto
      % splitData add_var function at generation time
      
      arg.name = name;
      arg.x    = x;
      arg.args = {varargin{:}};
      f.data(end+1) = arg;
    end
  end
end
