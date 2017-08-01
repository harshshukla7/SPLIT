function x = proj_negative(z)
%
% x = proj_positive(z)
%
% Project z onto the positive orthant
%

x = min(z,0);
