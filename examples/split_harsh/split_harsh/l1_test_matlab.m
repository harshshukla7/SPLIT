
    z=[-7, -5, -6];
    c=1;
    
  z = z/c;
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
x = (z+lam < 0) .* (z+lam)  + (z-lam > 0) .* (z-lam)
