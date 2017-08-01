classdef splitGen < handle
  % class for the generation of functions to do basic operations
  
  properties
    dat = []; % File to hold problem data
  end
  
  methods
    function gen = splitGen(dat)
      assert(isa(dat,'splitData'), 'Argument dat must be of type splitData')
      gen.dat = dat;
    end
    
    function mv_mult(gen, A, Astr, funcName, method)
      % Generate code for matrix-vector multiplication for a known matrix
      %
      % mv_mult(A, Astr, funcName, method)
      %
      % Creates c-code that implements the function y = A*x
      %   funcName(double y[size(A,1)], double x[size(A,2)]))
      %
      % method is one of: blas, dense, exhaustive, sparse
      % if unspecified, we guess as to the best approach
      
      density = nnz(A) / prod(size(A));
      sz      = max(size(A));
      
      if nargin < 5 % Guess the best method
        
        % Less than 100 rows, and 10% density, we do exhaustive
        if density < 0.1 && sz < 100
          method = 'exhaustive';
        elseif sz < 20
          method = 'exhaustive';
        elseif density < 0.3
          method = 'sparse';
        else
          method = 'blas';
        end
      end
      
      fprintf('Generating function %s using method %s\n', funcName, method);
      eval(sprintf('gen.gen_mv_%s(A, Astr, funcName);', method))
    end
    
    function solve(gen, A, Astr, funcName, method)
      % Generate code to solve a square system of linear equations A*x = b
      %
      % solve(A, Astr, funcName, method)
      %
      % Creates: funcName(double x[n], double b[n])
      %
      % Method = 'invert' or 'LDL'
      % Will guess a method if one isn't given
      
      gen.gen_solve_invert(A,Astr,[funcName '_invert']);

      gen.gen_solve_ldl(A, Astr, [funcName '_ldl']);

      
    end
  end
  
  methods (Hidden)
    function gen_solve_banded(gen, A, Astr, funcName)
      % Solve 
    end
    
    function gen_solve_banded_pd_symmetric(gen, A, Astr, funcName)
      % Solve the positive-definite symmetric system of equations by:
      % 1) Reducing bandwidth by permuting with symrcm
      % 2) Compute LU factorization (in c, with LAPACK)
      % Solve
      % a) Permute RHS (in c)
      % b) Solve (in c, with LAPACK)
      % c) Un-permute (in c)
      
%       % Compute permulation
%       p=symrcm(A);
%       gen.dat.add([Astr '_p'], p, 'int');
      
      % Produce a matrix of the banded format required by LAPACK
      % See dpbsv documentation for 'L' format
      %       Ab = A(p,p);
      Ab = A;
      KD = bandwidth(Ab); N = size(A,1);
      M  = zeros(KD,N);
      for i = 1:KD+1
        I = [i:N]; J = [1:N-i+1];
        dg = Ab(sub2ind(size(Ab),I,J));
        M(i,1:N-i+1) = dg;
      end      
      gen.dat.add(Astr, vec(M), 'double'); % Store column-wise
      gen.dat.define([Astr '_KD'], KD, 'int');
      
      pl('bool flag_%s_init = false;', funcName) % Used to indicate is the LU factorization has been done yet
      pl('void %s(double _x[%i], double _b[%i]) {', funcName, N, N)
      pl('char UPLO = ''L'';\nint  N    = %i;\nint  KD   = %i;\nint  NRHS = 1;\n',N,KD);
      pl('int  LDAB = %i+1;\nint  LDB  = %i;\nint  INFO;', KD, N);
      pl('  if (flag_%s_init == 0) {', funcName)
      pl('    // Factor matrix')
      pl('    dpbtrf_( &UPLO, &N, &KD, %s, &LDAB, &INFO );', Astr)
      pl('    flag_%s_init = true;', funcName)
      pl('  }')
      pl('  // Solve')
      pl('  memcpy(_x,_b,sizeof(double)*%i);', N)
      pl('  dpbtrs_( &UPLO, &N, &KD, &NRHS, %s, &LDAB, _x, &LDB, &INFO );', Astr)
      pl('}')      
    end
    
    function gen_solve_ldl(gen, A, Astr, funcName)
      % Solve the symmetric system of equations via LDL pre-factorization

      % Pre-factor the matrix A
      A = sparse(A);
      [L,D,p,S] = ldl(A,0.01,'vector');
      s = diag(S);
      s = s(p);

      % Create all the functions that we will need:
      
      % y  = L \ (s .* b(p));
      gen.lower_triangular_solve(L, [Astr '_L'], s, p, [funcName '_solve1']);

      % z  = inv(D) * y;
      gen.mv_mult(inv(D), 'iD', [funcName '_solve2'], 'exhaustive');

      % q  = L' \ z;
      % x(p,1)  = s .* q;
      gen.upper_triangular_solve(L', [Astr '_LT'], s, p, [funcName '_solve3']);
      
      % Create the function that ties everything together
      pl('void %s(double _x[%i], double _b[%i]) {', funcName, length(s), length(s))
      pl('  %s_solve1(_x,  _b);', funcName)
      pl('  %s_solve2(_b,  _x);', funcName)
      pl('  %s_solve3(_x,  _b);', funcName)
      pl('}')
    end
    
    
    function gen_solve_invert(gen, A, Astr, funcName)
      % Solve system of equations by pre-inverting A
      gen.mv_mult(inv(A), sprintf('%s_inv',Astr), funcName);
    end
    
    
    function gen_mv_blas(gen, A, Astr, funcName)
      % Multiply static matrix and vector using BLAS
      
      % Store data file
      gen.dat.add(Astr, vec(full(A)'));
      
      pl('void %s(double y[%i], const double x[%i]) {', funcName, size(A,1), size(A,2))
      pl('cblas_dgemv(CblasRowMajor, CblasNoTrans, %i, %i, 1.0, (double*)%s, %i, x, 1, 0.0, y, 1);',...
        size(A,1), size(A,2), Astr, size(A,1));
      
      pl('}')
    end
    
    function gen_mv_dense(gen, A, Astr, funcName)
      % Multiply static matrix and vector using for-loops
      
      % Store data file
      gen.dat.add(Astr, vec(full(A)'));
      
      pl('void %s(double y[%i], const double x[%i]) {', funcName, size(A,1), size(A,2))
      
      pl('for(int r=0; r<%i; r++) {', size(A,1))
      pl('  y[r] = 0.0;')
      pl('  for(int c=0; c<%i; c++) y[r] += %s[r*%i+c]*x[c];', size(A,2), Astr, size(A,1))
      pl('}')
      
      pl('}')
    end
    
    function gen_mv_exhaustive(~, A, ~, funcName)
      % Write out every non-zero multiplication
      %
      % Will be a stupid method for large matrices
      
      pl('void %s(double y[%i], const double x[%i]) {', funcName, size(A,1), size(A,2))
      
      A = full(A);
      
      for r=1:size(A,1)
        p('y[%i] = 0.0', r-1)
        
        if size(A,2) > 50
          [~,ind,val] = find(A(r,:));
          if ~isempty(ind)
            p('+ %.10g*x[%i]', [val;ind-1]);
          end
        else
          
          for c=1:size(A,2)
            if abs(A(r,c)) > 1e-6
              if A(r,c) == 1
                p(' + x[%i]', c-1)
              elseif A(r,c) == -1
                p(' - x[%i]', c-1)
              elseif A(r,c) < 0
                p(' - %.10g*x[%i]', abs(A(r,c)), c-1)
              elseif A(r,c) > 0
                p(' + %.10g*x[%i]', A(r,c), c-1)
              end
            end
          end
        end
        pl(';')
      end
      
      pl('}')
    end
    
    function gen_mv_sparse(gen, A, Astr, funcName)
      % Multiply static matrix stored in compressed form and vector
      
      % Write data in compressed format
      gen.dat.writeSparseMatrix(A,Astr);
      pl('void %s(double y[%i], const double x[%i]) {', funcName, size(A,1), size(A,2))
      
      pl('int i=0;')
      pl('for(int r=0; r<%i; r++) {', size(A,1))
      pl('  y[r] = 0.0;')
      pl('  for(int ic=0; ic<%s_nz_per_row[r]; ic++) {', Astr)
      pl('    y[r] += %s_dat[i]*x[%s_ind[i]];', Astr, Astr)
      pl('    i++;')
      pl('  }')
      pl('}')
      
      pl('}')
    end
  end
end
