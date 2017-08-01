function x = prox_quad(z, t, A, b, ~)
% 
% x = prox_quad(z, t, A, b, c)
%
% Compute the proximal operator of the quadratic function
%
% x = argmin 1/2*x'*A*x + b'*x + c + 1/(2t)||x-z||_2^2
%

n = size(A,1);
x = (eye(n) + t*A) \ (z-t*b);
