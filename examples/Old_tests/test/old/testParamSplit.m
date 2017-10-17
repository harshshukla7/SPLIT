clear
splitProb.clearProblem

x = splitvar(2);
y = parameter(2);

A = randn(10,2); b = ones(10,1); B = randn(10,2);
C = randn(8,2); D = randn(8,2);
c = ones(2,1);
e = randn(1,2); g = randn(1,2);

% Define a parametric problem
A*x <= b + B*y
% e*x == g*y + 1;
norm(C*x + D*y) <= c'*x + 5*sum(y);

minimize(x'*x + sum(y))

% Generate the problem data
prob = splitProb.genProblem

%%

% Solve the problem for a variety of different values of the parameter

for i = 1:1
 isfeasible = false;
 
 while ~isfeasible
  y.set(rand(2,1));
  
  % Confirm solution with YALMIP
  xx = sdpvar(2,1); yy = sdpvar(2,1);
  con = set(A*xx <= b + B*yy);% + set(e*xx == g*yy + 1);
  con = con + set(norm(C*xx + D*yy) <= c'*xx + 5*sum(yy));
  obj = xx'*xx + sum(yy);
  
  con = con + set(yy == y.val);
  
  ddd = solvesdp(con, obj);
  if ddd.problem == 0
   isfeasible = 1;
  end
 end
 
 [sol, stats, data_to_export] = admm(prob);
 splitProb.setSolution(prob, sol); 
 
 fprintf('Error = %e\n', norm(x.val - double(xx)));
 
end

return

%%

clear

% Must always call clearProblem before starting!
% (or use *clear all*, clear is not enough)
splitProb.clearProblem;

N = 5;
n = 3;
p = 2;
m = 2;
[A,B,C,D] = ssdata(drss(n,p,m));

x = splitvar(n, N);
u = splitvar(m, N-1);
P = splitvar(n, n, 'symmetric');
Q = randn(n); Q = Q*Q';

x0 = parameter(n,1);
x(:,1) = x0;

for i = 1:N-1
 x(:,i+1) == A*x(:,i) + B*u(:,i);
end

% Box constraints
u.inBox(-ones(m,N-1), ones(m,N-1));

% Upper / lower bounds
-5 <= x <= 5;

% Cost must be linear for now
minimize(ones(n,1)'*x(:,end));

%% Generate the problem data
prob = splitProb.genParametricProblem;

%% Solve the problem for various initial states

X0 = randn(n,10);
for k = 1:size(X0,2)
 set(x0, X0(:,k));
 [sol, stats] = admm(prob);
 
 if sol.problem
  error('Could not solve the problem with ADMM - infeasible?')
 end
 
 % Recover the solution
 splitProb.setSolution(prob, sol);
 
 % Solution is now stored in the variables
 cprintf('_blue', '                OPTIMAL SOLUTION               \n')
 fprintf('x = \n'); disp(full(x.val))
 fprintf('u = \n'); disp(full(u.val))
 
 % Compare to Yalmip
 
 yal_x = sdpvar(n, N);
 yal_u = sdpvar(m, N-1);
 
 con = [];
 for i = 1:N-1
  con = con + set(yal_x(:,i+1) == A*yal_x(:,i) + B*yal_u(:,i));
 end
 
 con = con + set( -5 <= yal_x <= 5) + set(-1 <= yal_u <= 1);
 con = con + set(yal_x(:,1) == X0(:,k));
 
 obj = ones(n,1)'*yal_x(:,end);
 
 result = solvesdp(con, obj, sdpsettings('verbose', 0));
 
 if result.problem, error('YALMIP could not solve problem'); end
 
 % Solution is now stored in the variables
 cprintf('_blue', '                ERROR TO YALMIP               \n')
 fprintf('%15s = %.2e\n', '||x - yal_x||', norm(x.val-double(yal_x)));
 fprintf('%15s = %.2e\n', '||u - yal_u||', norm(u.val-double(yal_u)));
end


