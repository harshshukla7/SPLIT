%--------------------------------------------------
% blkdiag
%--------------------------------------------------
function nvar = blkdiag(varargin)
nvar = varargin{1};
for i = 2:length(varargin)
  t = varargin{i};
  Z1 = sparse(size(nvar,1), size(t,2));
  Z2 = sparse(size(t,1), size(nvar,2));
  nvar = [nvar Z1; Z2 t];
end
