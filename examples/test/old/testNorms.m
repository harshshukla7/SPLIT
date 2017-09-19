clear all

%%

Q = randn(2); Q = Q*Q';

x = splitvar(2);
f = norm(x,1);

% x'*Q*x <= 3;
% 2*f^2 + x'*Q*x <= 4.3;
% x'*Q/3*x <= 3.1;
% 
% f^2 + norm(x,1) <= 3;

% 0.3*norm(randn(2)*x,1) + x'*x + sum(abs(x)) + max(x)^2/2 <= sum(randn(2)*x + 3);
% 0.3*norm(randn(2)*x,1) + x'*x + sum(x) + x(2)^2/2 <= sum(randn(2)*x + 3);

norm(Q*x + 3,inf) < 3 + sum(x);
x(1)^2 + x(2) <= 5;

% X = splitvar(3,3,'symmetric');

% f1 = norm(x,1);
% f2 = norm(x,2);
% finf = norm(x,inf);
% finf2 = norm(x,'inf');

% ab = abs(x);
% z = splitvar(2,3);
% mx = max(x,2*x);
% mn = min(randn(size(x,1),1), z(:,2));

% norm(randn(4,2)*x + randn(4,1), inf) <= 3;
% 7*f1^2 + 3*f2 + 3.4*mx + 6*x(1)^2 <= 4;
% finf/4 < 5;
% 4*f2^2/2 <= 6;
% 3*x'*Q*x <= [4;5]'*x + 3;
% 3.4*norm(randn(3,2)*x + ones(3,1)) <= [5;6]'*x + 4;
% -1 <= x <= 3;
% X - f1 >= 0;

% y = splitvar(2,5);
% u = splitvar(1,4);
% A = randn(2);
% B = randn(2,1);
% for i = 1:4
%   y(:,i+1) == A*y(:,i) + B*u(:,i);
% end
% 
% 
% t = randn(1,2)*x + 7*f1^2 + 3*f2 + 3.4*x(2) + 6*x(1)^2 + randn(1,2)*z(:,2);
% sp = splitProb.instance;
% 
% t = [[randn(1,2)*x + 7*f1^2 + 3*f2 + 3.4*x(2) + 6*x(1)^2 + randn(1,2)*z(:,2);randn(1,2)*x] ...
%  [x(1);min(mn + x(2))]];

return

y = splitvar(2,5);
u = splitvar(1,4);
A = randn(2);
B = randn(2,1);
for i = 1:4
  y(:,i+1) == A*y(:,i) + B*u(:,i);
end




sp = splitProb.instance;
con=sp.constraints(1);
a = con.a(:); b = con.b(:);
i=1;

return


%%

%
% Test the handling of norm functions

x = splitvar(2);

% b = randn(3,1);
% norm(x, 1) + norm(x, 'inf') + 2*norm(x) <= 3;
% a = randn(2,1);
% a'*x == 0;

max(x(1), x(1)-x(2)) <= 3;
x(2) >= 0.1;
t = sum(abs(sum(x,2)));
t <= -norm(x(1),1) + 3;

obj = -sum(x) + 0.1*norm(x,1);
minimize(obj);


prob = splitProb.genProblem;

%%

[sol, stats] = admm(prob)

splitProb.setSolution(prob, sol);

%%

% Compare to YALMIP

y = sdpvar(2,1);
con = set(max(y(1), y(1)-y(2)) <= 3);
con = con + set(y(2) >= 0.1);
q = sum(abs(sum(y,2)));
con = con + set(q <= -norm(y(1),1) + 3);

yObj = -sum(y) + 0.1*norm(y,1);
solvesdp(con, yObj, sdpsettings('solver','linprog'))

fprintf('-----------------------------------\n');
fprintf('Comparison to YALMIP \n');
fprintf('  Optimal value error = %.2e\n', obj.val - double(yObj));
fprintf('  optimizer error     = %.2e\n', norm(x.val - double(y)))
fprintf('-----------------------------------\n');
