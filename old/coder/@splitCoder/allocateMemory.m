function allocateMemory(cdr)
%
% Allocate memory for internal variables
%
% Must be called before any of the internal variables are used!
%

if ~isempty(cdr.allocVars)
  % Declare all internal variables as persistent
  cdr.print('persistent ');
  for i = 1:length(cdr.allocVars)
    cdr.print([cdr.allocVars(i).name ' '])
  end
  cdr.print('\n')
  
  % Allocate memory
  cdr.print('if isempty(%s)\n', cdr.allocVars(1).name)
  for i = 1:length(cdr.allocVars)
    if ischar(cdr.allocVars(i).m), m = cdr.allocVars(i).m;
    else                           m = sprintf('%i', cdr.allocVars(i).m);
    end
    if ischar(cdr.allocVars(i).n), n = cdr.allocVars(i).n;
    else                           n = sprintf('%i', cdr.allocVars(i).n);
    end
    
    cdr.print('\t%s = zeros(%s, %s);\n', cdr.allocVars(i).name, m, n)
  end
  cdr.print('end\n')
end

