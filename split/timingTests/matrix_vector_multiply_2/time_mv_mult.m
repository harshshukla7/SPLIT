clear all

% Generating timings for the multiplication of a matrix and a vector as a
% function of sparsity and size

% Comparison of computational techniques:
% 1. Exhaustive code-generation
% 2. Dense matrix-vec multiplication
% 3. BLAS
% 4. Sparse matrix storage

global f

methods = {'exhaustive','dense','sparse','blas'};

fDat = fopen('data.dat','w+');
fprintf(fDat, 'Name, Rows, Density, Time\n');
fclose(fDat);

KK = ceil(logspace(1,3,30));
DD = logspace(-2,0,5); DD = ceil(DD*100)/100;
testNum = 1;
for d = DD
  for k = KK
    n = k;
    A = sprandn(n,n,d);
    x = randn(n,1);
    sol = A*x;
    
    fprintf('----- n = %i, d = %f, nnz = %i ----- (%i of %i)\n', n, d, nnz(A), testNum, length(KK)*length(DD));
    testNum = testNum + 1;
    
    dat = splitData('splitData');
    gen = Gen(dat);
    
    f = dat.fHdr;
    
    dat.define('n', n);
    dat.define('density', d);
    
    dat.add('x',x);
    dat.add('sol',sol);

    if n < 100
      gen.mv_mult(A, 'A_exhaustive', 'mv_exhaustive', 'exhaustive');
    else
      dat.define('mv_exhaustive(y,x)','');
    end
    gen.mv_mult(A, 'A_sparse', 'mv_sparse', 'sparse');
    gen.mv_mult(A, 'A_dense', 'mv_dense', 'dense');
    gen.mv_mult(A, 'A_blas', 'mv_blas', 'blas');
    
    dat.writeFile;
    
    [q,out] = system('gcc -O3 -framework Accelerate time_mv_mult.c splitTimer.c splitLoad.c');
    if q ~= 0, error('Could not compile!'); end
    system('./a.out');
  end
end
