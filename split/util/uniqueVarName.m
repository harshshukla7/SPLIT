function name = uniqueVarName
% name = uniqueVarName
%
% Generates a unique variable name
%

persistent id
if isempty(id)
  id = 1;
end

name = sprintf('_var_%i_', id);
id = id + 1;
