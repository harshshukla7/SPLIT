%--------------------------------------------------
% Horizontal Concantenation
%--------------------------------------------------
function nvar = horzcat(varargin)

% Drop any empty cells
varargin = {varargin{~cellfun('isempty',varargin)}};
varargin = {varargin{~cellfun(@(x) length(x(:))==0, varargin)}};

n = size(varargin{1},1);
basis = sparse(splitProb.numVars+1,0);
for i = 1:length(varargin)
  t = varargin{i};
  
  assert(n == size(t, 1), 'splitvar:sizechk', 'Matrix dimensions must agree')
  if isnumeric(t), t = splitvar(size(t), [], t(:)'); end
  b = t.basis;
  basis(:,end+1:end+size(b,2)) = b;
end
nvar = splitvar(n, size(basis,2)/n, basis);
