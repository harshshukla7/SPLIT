function mldivide(cdr, K, x, y, mode)
%
% cdr.mldivide(K, x, y, mode)
%
% Solve the linear system:
%   y = K\x
%
% K : Constant matrix
% x : [string] 
% y : [string] Adds a working variable called y
%
% mode is one of:
%  - LU           Generate explicit LU decomposition
%  - Invert KKT   Pre-compute and store the inverse of K

cdr.addVar(y, size(K,1), 1);

assert(size(K,1) == size(K,2), 'Matrix K must be square')
assert(rank(K) == size(K,1), 'K must be full rank')

switch lower(mode)
  case 'lu'
    error('Have not implemented explicit LU decomposition yet')
  case 'invert kkt'
    iName = sprintf('iK_%i', randi(1000));
    cdr.addConst(iName,'Pre-computed matrix inverse')
    eval(sprintf('cdr.%s = inv(K);', iName));
    cdr.print('%s = %s*%s;\n', y, iName, x)
  otherwise
    error('Unknown mldivide solve mode %s', mode)
end


% %% >>>>>>>>>>> RUN >>>>>>>>>>>
% if strcmp(cdr.KKTSolveMode, 'LU')
%   cdr.print('function x = solveKKT(b, x)\n');
%   
%   % Allocate memory for the LU temporary variable
%   cdr.print('persistent y\n');
%   cdr.print('if isempty(y), y = zeros(%i,1); end\n', size(K,1));
%   
%   % Compute the LU factorization
%   [LL UU p q] = lu(sparse(K),'vector');
%   cdr.genSparseLU(LL,UU,p,q,'x','y','b');
%   
%   cdr.print('end\n');
% end
% 
