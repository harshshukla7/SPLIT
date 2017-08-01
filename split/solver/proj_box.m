function x = proj_box(z, lb, ub)
%
% x = proj_box(z, lb, ub)
%
%  x = argmin ||x - z|| s.t. lb <= z <= ub
%

x = min([max([z lb],[],2) ub],[],2);
