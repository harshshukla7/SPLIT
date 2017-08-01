function x = proj_quadball(z, c)
%
% x = proj_quadball(z)
%
% Project z onto the quadratic ball
% 
% Implementation details: Colin's code proj_normBall : case 2

nz = norm(z);

% We do not need to compute square root of c because in
% the class ellipseConSet.m is already computed

if nz <= c
   x = z;
else
   x = z * c / nz;
end
