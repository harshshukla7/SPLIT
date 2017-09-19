function x = prox_norm(z, t, p)
%
% x = prox_norm(z, t, p)
%
% Compute the proximal operator for the p-norm
%
% x = argmin t*||x||_p + 1/2*||z - x||_2^2 
%


% prox_tf(v) = v - t*proj(v/t)

switch p
  case 1
    pDual = inf;
  case 2
    pDual = 2;
  case inf
    pDual = 1;
  otherwise
    error('Cannot handle the %s-norm', p)
end

x = z - t*proj_normBall(z/t, pDual);
