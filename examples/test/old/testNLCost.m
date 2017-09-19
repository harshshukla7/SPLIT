clear all

fprintf('---> Solving with split\n');
x = splitvar(2)

% Test processing of non-linear objective functions
obj = 3*x'*x + 1/3.4*norm(x) + 3*sum(x) + [6 7]*x;

1 <= x <= 4;

minimize(obj)

prob = splitProb.genProblem

sol = admm(prob)

splitProb.setSolution(prob, sol)

%%
fprintf('---> Solving with yalmip\n');
y = sdpvar(2,1);
solvesdp(set(1 <= y <= 4), 3*y'*y + 1/3.4*norm(y) + 3*sum(y) + [6 7]*y)

fprintf('Error = %.4e\n', norm(x.val - double(y)))
