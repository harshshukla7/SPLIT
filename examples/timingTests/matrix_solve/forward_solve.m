function funcname = forward_solve(cdr, L, Lstr, p)
% funcname = forward_solve(cdr, L, Lstr, p)
%
% Generate sparse code for forward solve of 
%   L*x = P'*S*b
%
% L    : lower triangular matrix
% Lstr : name of the matrix L (for dense generation)
% P    : permutation matrix
% S    : diagonal scaling matrix
%
% Returns the name of the generated function
%
% funcname(double *x, double *b)
%

funcname = sprintf('forward_solve_%i', cdr.getID());

cdr.pl;
cdr.pl('void %s(double x[%i], double b[%i])', funcname, size(L,2), size(L,1))
cdr.pl('{')

% if density > cdr.sparseProdLimit
  % Produce dense code
  
  % Store the strictly lower triangular matrix row-wise
  lvec = [];
  for i = 1:size(L,1)
    lvec = [lvec L(i,1:i-1)];
  end
  % Store the inverse of the diagonal
  iDiagL = 1./diag(L);
  
  LTristr  = cdr.define_constant(sprintf('%s_ltri',Lstr), lvec);
  LDiagstr = cdr.define_constant(sprintf('%s_diag',Lstr), lvec);

  % Generate forward solve code
  cdr.pl('int i=0;')
  cdr.pl('for(int r=0; r<%i; r++) {', size(L,1))
  cdr.pl('  x[r] = b[r];')
  cdr.pl('  for(int c=0; c<r; c++) x[r] -= %s[i++]*b[c];', LTristr)
  cdr.pl('  x[r] *= %s[r];', LDiagstr);
  cdr.pl('}')
  
% else
% end

cdr.pl('}')
