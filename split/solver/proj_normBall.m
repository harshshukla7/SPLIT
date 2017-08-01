function x = proj_normBall(z, p, c)
%
% x = proj_normBall(z, p, c)
%
% Compute the projection onto the p-norm-ball
%
% x = argmin ||z - x||_2 s.t. ||x||_p <= c
%
% if c is not specified, then it's assumed to be 1
%


% Projections from Boyd's lecture notes
% http://www.seas.ucla.edu/~vandenbe/236C/lectures/proxop.pdf

if nargin < 3, c = 1; end

switch p
  case 1
    if c ~= 1, z = z/c; end
    
    nz = norm(z, 1);
    if nz <= 1
      lam = 0;
    else
      % lam solves the equation sum max(|x|-lam, 0) == 1
      v = sort(abs(z));
      for i = 1:length(v)
        if sum(max(v-v(i),0)) < 1, break; end
      end
      % We know that lam lies between v(i-1) and v(i), and that i > 1
      % Solve for lam
      lam = (sum(v(i:end))-1)/(length(v)-i+1);
    end
    
    x = (z+lam < 0) .* (z+lam)  + (z-lam > 0) .* (z-lam);
    
  case 2
      
    nz = norm(z);
    if nz <= c
      x = z;
    else
      x = z * c / nz;
    end
    
  case inf
    x = proj_box(z, -c*ones(size(z)), c*ones(size(z)));
    
  otherwise
    error('Unknown norm-ball %s', p)
end