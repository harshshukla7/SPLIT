function addVar(cdr, name, m, n)
%
%  cdr.addVar(name, m, n)
%
% Add an internal variable
%
% Causes coder to allocate memory for it
%
% name = variable name
% m,n  = size. Can be strings or constants
%

% Check if we've already got this variable
if ~isempty(cdr.allocVars)
  match = strcmp({cdr.allocVars.name}, name);
else
  match = false;
end

if any(match)
  i = find(match, 1);
  var = cdr.allocVars(i);
  
  % Confirm that the sizes match
  assert(ischar(m) == ischar(var.m), 'Inconsistent re-definition of internal variable %s', name)
  if ischar(m)
    assert(strcmp(m, var.m), 'Inconsistent re-definition of internal variable %s', name)
  else
    assert(m == var.m, 'Inconsistent re-definition of internal variable %s', name)
  end
  
  assert(ischar(n) == ischar(var.n), 'Inconsistent re-definition of internal variable %s', name)
  if ischar(n)
    assert(strcmp(n, var.n), 'Inconsistent re-definition of internal variable %s', name)
  else
    assert(n == var.n, 'Inconsistent re-definition of internal variable %s', name)
  end
else
  len = length(cdr.allocVars) + 1;
  
  cdr.allocVars(len).name = name;
  cdr.allocVars(len).m    = m;
  cdr.allocVars(len).n    = n;
end

