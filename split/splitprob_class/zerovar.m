function var = zerovar(n,m)
%
% Create a split variable that is all zero
%

if nargin < 2, m = 1; end
if nargin < 1, n = 1; end

if numel(n) == 2, sz = n;
else              sz = [n m];
end

var = splitvar(n,m,spalloc(1,prod(sz),0));
