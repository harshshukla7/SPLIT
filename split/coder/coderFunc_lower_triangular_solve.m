classdef coderFunc_lower_triangular_solve < coderFunc
  
  methods
    function f = coderFunc_lower_triangular_solve(L, Lstr, s, p, funcname)
      % f = coderFunc_lower_triangular_solve(L, Lstr, s, p, funcname)
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

      f = f@coderFunc('void %s(real _x[%i], real _b[%i])', funcname, size(L,2), size(L,1));
      
      % Store the inverse of the diagonal
      iLDiag = sprintf('%s_iDiag', Lstr);
      f.add_var(iLDiag, 1./diag(L));
      
      % Store the lower-triangular portion of L in sparse form
      Ltri = sprintf('%s_lower_tri', Lstr);
      f.add_var(Ltri, sparse(tril(L,-1)));
      
      % % Store the permutation vector and scaling
      pName = sprintf('%s_p', Lstr);
      sName = sprintf('%s_s', Lstr);
      f.add_var(pName, p-1, 'type', 'int', 'storage_method', 'dense');
      f.add_var(sName, s,   'type', 'real', 'storage_method', 'dense');
      
      f.pl('int i = 0;')
      f.pl('for(int r=0; r<%i; r++) {', size(L,1));
      f.pl('  _x[r] = %s[r] * _b[%s[r]];', sName, pName);
      f.pl('  for(int j=0; j<%s_nz_per_row[r]; j++) {', Ltri);
      f.pl('    _x[r] -= %s_dat[i] * _x[%s_ind[i]];', Ltri, Ltri);
      f.pl('    i++;');
      f.pl('  }');
      f.pl('  _x[r] *= %s[r];', iLDiag);
      f.pl('}');
    end
  end
end
