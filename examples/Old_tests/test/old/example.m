clear all

x = splitvar(2);
f = norm(x,1);

x'*Q*x <= 3;
2*f^2 + x'*Q*x <= 4.3;
f^2 + norm(x,1) <= 3;

norm(Q*x + 3,inf) < 3 + sum(x);
x(1)^2 + x(2) <= 5;




% Optimization variables
x = splitvar(n, N);
u = splitvar(m, N-1);
x(:,1) = parameter(n,1);

% Objective and dynamics
obj = 0;
for i = 1:N-1
  x(:,i+1) == A*x(:,i) + B*u(:,i);
  obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i);
end
obj = obj + x(:,end)'*x(:,end);

% set up constraints to generate the problem 
-5 <= x <= 5;
-1 <= u <= 1;
norm(x(2,:), 2) + x(:)'*x(:) <= 4;

minimize(obj);

% Generate problem
prob = splitProb.genProblem;

% Solve in matlab
[solution, stats] = split_solve(prob);

% Generate c-code
splitCodegen('code.c');
