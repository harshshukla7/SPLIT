function z = proj_positive(z,~)
%
% z = proj_positive(z)
%
% Project z onto the positive orthant
%

% x = max(z,0);
% z(z < 0) = 0;
for(i=1:length(z))if(z(i)<0)z(i)=0;end;end
