%--------------------------------------------------
% Vertical concantenation
%--------------------------------------------------
function nvar = vertcat(nvar, varargin)
nvar = nvar';
for i = 1:length(varargin)
  nvar = [nvar varargin{i}'];
end
nvar = nvar';
