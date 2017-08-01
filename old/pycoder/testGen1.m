clear all

setenv('PATH',['/Users/cnjones/anaconda/bin:', getenv('PATH')])

A = [1 1;1 0];
B = [1;0.5];
N = 5;

x = splitvar(2,N);
u = splitvar(1,N-1);
x0 = parameter(2,1);

x(:,2:end) == A*x(:,1:end-1) + B*u;
x(:,1) == x0;

-5 <= x(:) <= 5;
-1 <= u(:) <= 1;

% Add something more complex
norm(x(:,3),2) + norm(x(:,4),inf) + norm(x(:,5),1) <= 4;

minimize(x(:)'*x(:) + u(:)'*u(:))

prob = splitProb.genProblem;
cdr  = pyCoderData(prob);
cdr.genData

!./runBuild


mex admm_mex.c admm.c matrix_ops.c splitTimer.c

%% Now check the solution

for i = 1:10
  par = 2*randn(2,1);
  x0.set(par);
  
  fprintf('===> Solving in matlab\n');
  sol = admm(prob);
  
  fprintf('===> Solving in c via mex\n');
  y = admm_mex(par);
  
  [x,~,~,itr,tm] = funcADMM_mex(par);
  fprintf('===> Solving in c via mex : %e\n', tm/itr);
  
  
  err = norm(y-sol.x);
  fprintf('Norm difference between c solution and m solution: %e\n', err);
  
end
