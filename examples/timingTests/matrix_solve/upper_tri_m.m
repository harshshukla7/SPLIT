function x = upper_tri_m(U, s, p, z)

close all
dat = splitData('blah');
[len,I,vec] = writeSparseMatrix(dat, U, 'U');

i = lenght(I);
for r = size(U,1) : -1 : 1
  x(r) = z(r);
  
  for j=r+1 : len(r)
  
  for(int j=r; j<%s_nz_per_row[r]; j++) {', Utri);
pl('      _x[r] -= %s_dat[i] * _x[%s_ind[i]];', Utri, Utri);
pl('      i--;');

%   for j = r+1 : size(U,2)
%     x(r) = x(r) - U(r,j)*x(j);
%   end
end
