clear all
close all

%%
% Produce an example matrix to solve from a random MPC problem

n = 10;
p = 1;
m = 2;
[A,B,C,D] = ssdata(drss(n,p,m));

NN = 2:2:50;
% NN = 50;

fDat = fopen('data.dat','w+');
fprintf(fDat, 'Name, Rows, Density, Horizon, Time\n');
fclose(fDat);

for N = NN
  fprintf('Horizon %i of %i\n', N, max(NN));
  
  % Setup a simple MPC problem
  % Upper/lower bounds and an ellipsoidal terminal set
  splitProb.clearProblem;
  x = splitvar(n,N);
  u = splitvar(m,N-1);
  x(:,2:end) == A*x(:,1:end-1) + B*u;
  -10 <= x(:) <= 10;
  -1 <= u(:) <= 1;
  P = randn(n); P = P*P';
  x(:,end)'*P*x(:,end) <= 1;
  x(:,1) == randn(n,1);
  
  minimize(x(:)'*x(:) + 0.1*u(:)'*u(:));
  
  prob = splitProb.genProblem;
  [sol, stats_out] = admm(prob);
  KKT = stats_out.KKT;
  
  %%%% For now - test symmetric form
  KKT = KKT*KKT';
  
  % Solve
  dd = nnz(KKT) / prod(size(KKT));
  
  dat = splitData('splitData');
  global f
  f = dat.fHdr;
  gen = splitGen(dat);
  
  gen.solve(KKT, 'A', 'solveKKT');

  p=symrcm(KKT);
  gen.gen_solve_banded(KKT(p,p), 'A_bnd', 'solveKKT_banded');
  
  b = randn(size(KKT,1),1);
  dat.add('b', b);
  
  dat.add('sol', KKT\b);
  dat.define('n',size(KKT,1),'int');
  dat.define('density',dd,'float');
  dat.define('horizon',N,'int');
  
  dat.writeFile;
  
  [q,out] = system('gcc -O3 -framework Accelerate test_matrix_solve.c splitTimer.c splitLoad.c');
  if q ~= 0, error('Could not compile!'); end
  system('./a.out');
end
return


%%

KK = ceil(logspace(1,3,30));
DD = logspace(-3,0,5); DD = ceil(DD*100)/100;
testNum = 1;
for d = DD
  for k = KK
    fprintf('Density = %f, size = %i\n', d, k);
    
    n = k;
    
    A = [];
    while rank(full(A)) < n
      A = sprandn(n,n,d);
      A = triu(A);
      A = A + A';
      
      A = A + diag(rand(n,1))*1e-3;
    end
    dd = nnz(A) / prod(size(A));
    
    dat = splitData('splitData');
    global f
    f = dat.fHdr;
    gen = splitGen(dat);
    
    gen.solve(A, 'A', 'solveKKT');
    
    b = randn(n,1);
    dat.add('b', b);
    
    dat.add('sol', A\b);
    dat.define('n',n,'int');
    dat.define('density',d,'float');
    
    dat.writeFile;
    
    [q,out] = system('gcc -O3 -framework Accelerate test_matrix_solve.c splitTimer.c splitLoad.c');
    if q ~= 0, error('Could not compile!'); end
    system('./a.out');
  end
end


% [L,D,p,S] = ldl(A,0.01,'vector');
% s = diag(S);
% s = s(p);
%
% x  = L \ (s .* b(p));
% y  = D \ x;
% q  = L' \ y;
% z(p,1)  = s .* q;
%
% sol = inv(A)*b;
%

% full(y)'
