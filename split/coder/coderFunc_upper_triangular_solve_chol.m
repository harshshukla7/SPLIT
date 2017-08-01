classdef coderFunc_upper_triangular_solve_chol < coderFunc
  
  methods
    function f = coderFunc_upper_triangular_solve_chol(U, Ustr, p, funcname)
      % funcname = upper_triangular_solve(U, Ustr, s, p, funcname)
      %
      % Generate sparse code for the solve of
      %   x(p) = (inv(U) * z)
      %
      % U    : upper triangular matrix
      % Ustr : name of the matrix U
      % p    : permutation vector
     
      %
      % funcname(double *x, double *z)
      %
      % NOTE : z is modified after call
      
      %       q  = L' \ z;
      %       x(p,1)  =  q;
      
      f = f@coderFunc('void %s(double _x[%i], double _z[%i])', funcname, size(U,2), size(U,1));
      
      % Store the inverse of the diagonal
      iUDiag = sprintf('%s_iDiag', Ustr);
      f.add_var(iUDiag, 1./diag(U));
      
      % Store the upper-triangular portion of U in sparse form
      Utri = sprintf('%s_upper_tri', Ustr);
      f.add_var(Utri, sparse(triu(U,1)));
      
      % Store the scaling vector
      %sName = sprintf('%s_s', Ustr);
      %f.add_var(sName, full(s));
      
      f.pl('int i = %s_ind_len-1;', Utri)
      f.pl('for(int r=%i; r>=0; r--) {', size(U,1)-1);
      f.pl('  _x[r] = _z[r];');
      f.pl('  for(int j=0; j<%s_nz_per_row[r]; j++) {', Utri);
      f.pl('    _x[r] -= %s_dat[i] * _x[%s_ind[i]];', Utri, Utri);
      f.pl('    i--;');
      f.pl('  }');
      f.pl('  _x[r] *= %s[r];', iUDiag);
      f.pl('}');
      
      % Scale and permute
      for i = 1:size(U,1)
        f.pl('  _z[%i] =  _x[%i];', p(i)-1, i-1);
      end
      
      % Copy to output
      f.pl('  memcpy(_x, _z, sizeof(double)*%i);', size(U,1));
    end
  end
end
