function lower_triangular_solve(gen, L, Lstr, s, p, funcname)
% funcname = forward_solve(gen, L, Lstr, s, p)
%
% Generate sparse code for forward solve of 
%   L*x = s .* b(p)
%
% L    : lower triangular matrix
% Lstr : name of the matrix L (for dense generation)
% p    : permutation vector
% s    : scaling vector
%
% funcname(double *x, double *b)
%
% NOTE : b is changed on exit!
%

dat = gen.dat;

% Store the inverse of the diagonal
iLDiag = sprintf('%s_iDiag', Lstr);
dat.add(iLDiag, 1./diag(L));

% Store the lower-triangular portion of L in sparse form
Ltri = sprintf('%s_lower_tri', Lstr);
dat.writeSparseMatrix(tril(L,-1), Ltri);

% % Store the permutation vector and scaling
pName = sprintf('%s_p', Lstr);
sName = sprintf('%s_s', Lstr);
dat.add(pName, p-1, 'unsigned int');
dat.add(sName, s);

pl;
pl('void %s(double _x[%i], double _b[%i])', funcname, size(L,2), size(L,1))
pl('{')
pl('int i = 0;')
pl('for(int r=0; r<%i; r++) {', size(L,1));
pl('  _x[r] = %s[r] * _b[%s[r]];', sName, pName);
pl('  for(int j=0; j<%s_nz_per_row[r]; j++) {', Ltri);
pl('    _x[r] -= %s_dat[i] * _x[%s_ind[i]];', Ltri, Ltri);
pl('    i++;');
pl('  }');
pl('  _x[r] *= %s[r];', iLDiag);
pl('}');
pl('}');
pl;
