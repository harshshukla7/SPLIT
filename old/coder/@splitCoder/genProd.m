function genProd(cdr,y,A,sA,x)
%
% cdr.genProd(y,A,sA,x)
%
% Generates code for y = A*x
%
% y = string (name of output)
% A = matrix A (data)
% A = string (name of matrix A)
% x = string (name of input)
%
% Creates an internal variable called 'y'

if nnz(A) / numel(A) <= cdr.MatVec_SparseLimit
  % Generate code for the sparse function
  funcName = sprintf('%s_func_%i',y,randi(1000));

  % Function to call now
  cdr.print('%s = %s(%s, %s);\n', y, funcName, x, y);
  
  % Queue this function to be generated later
  cdr.addFunction(sprintf('function y = %s(x, y)', funcName), ...
    cdr.genSparseMatVec(A, 'x', 'y'),...
    sprintf('Multiplying matrix %s by vector %s', sA, x));
else
  % Just compute the product directly here
  cdr.print('%s = %s * %s;\n', y, sA, x)
end

% Add an internal variable to store the result
cdr.addVar(y, size(A,1), 1);
