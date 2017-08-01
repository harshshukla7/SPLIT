clear all

N = 5;
% A = [2 -1;1 0];
% B = [1;0];

n = 4;
m = 2;

[A,B,C,D] = ssdata(drss(n,m,m));

x  = splitvar(n,N);
u  = splitvar(m,N-1);
x0 = parameter(n,1);

x(:,2:end) = A*x(:,1:end-1) + B*u;
x(:,1) == x0;

-5 <= x(:) <= 5;
-1 <= u(:) <= 1;
 
for i = 1:N
  norm(x(:,i)) + norm(x(:,i),inf) <= 4;
end

minimize(x(:)'*x(:) + u(:)'*u(:));

prob = splitProb.genProblem;
cdr  = splitcoder(prob);

cdr.gen;
