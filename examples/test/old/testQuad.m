clear

splitProb.clearProblem

x = splitvar(3,1);
y = splitvar(2,1);

% Form a problem with various quadratic terms

Q = randn(3); Q= Q*Q';
H = randn(2); H=-H*H';

norm(x, 1) + 2*norm(x+randn(3,2)*y, inf) <= 5;
% y'*y <= 2;

% minimize(x'*Q*x - y'*H*y);
minimize(norm(Q*x,1))

prob = splitProb.flatten;

%% Generate the problem data
% prob = splitProb.genProblem;

% %% Solve the problem
% [sol,stats] = admm(prob);
% 
% if sol.problem
%   error('Could not solve the problem with ADMM - infeasible?')
% end
% 
% clf
% semilogy([stats.rDual], 'b');
% hold on; grid on
% semilogy([stats.rPrimal], 'm')
% legend('Dual residual', 'Primal residual')
