function upper_triangular_solve(gen, U, Ustr, s, p, funcname)
% funcname = upper_triangular_solve(gen, U, Ustr, s, p, funcname)
%
% Generate sparse code for the solve of 
%   x(p) = s .* (inv(U) * z)
%
% U    : upper triangular matrix
% Ustr : name of the matrix U
% p    : permutation vector
% s    : scaling vector
%
% funcname(double *x, double *z)
%
% NOTE : z is modified after call

%       q  = L' \ z;
%       x(p,1)  = s .* q;

dat = gen.dat;

% Store the inverse of the diagonal
iUDiag = sprintf('%s_iDiag', Ustr);
dat.add(iUDiag, 1./diag(U));

% Store the upper-triangular portion of U in sparse form
Utri = sprintf('%s_upper_tri', Ustr);
dat.writeSparseMatrix(triu(U,1), Utri);

% % Store the scaling vector
sName = sprintf('%s_s', Ustr);
dat.add(sName, s);

pl;
pl('void %s(double _x[%i], double _z[%i])', funcname, size(U,2), size(U,1))
pl('{')
pl('  int i = %s_ind_len-1;', Utri)
pl('  for(int r=%i; r>=0; r--) {', size(U,1)-1);
pl('    _x[r] = _z[r];');
pl('    for(int j=0; j<%s_nz_per_row[r]; j++) {', Utri);
pl('      _x[r] -= %s_dat[i] * _x[%s_ind[i]];', Utri, Utri);
pl('      i--;');
pl('    }');
pl('  }');

% Scale and permute
for i = 1:size(U,1)
pl('  _z[%i] = %s[%i] * _x[%i];', p(i)-1, sName, i-1, i-1);
end

% Copy to output
pl('  memcpy(_x, _z, sizeof(double)*%i);', size(U,1));

pl('}');
pl;
