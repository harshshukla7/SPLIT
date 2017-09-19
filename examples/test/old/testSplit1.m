clear all
%%

% Test simple operations in split

x = splitvar(2,1);

norm(x - 2*ones(2,1),inf) <= 1;

minimize(x'*x)

%%

prob = splitProb.genProblem;

%%
sol = admm(prob)
