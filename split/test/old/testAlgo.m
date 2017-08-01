clear 

algo = splitAlgo;

load ADMM_test_data

x   = algoVar(n,1, 'x');
y   = algoVar(sum(nProxVars),1,'y');
lam = algoVar(sum(nProxVars),1,'lam');
q   = algoVar(size(L,1),1,'q');

tmp = algoVar(size(KKT,1), 1, 'tmp');

rho = algoVar(1,1,'rho');

% while ...

algo.comment('Test: tuning the step-size rho: increase rho by 0.99');
algo.assign(rho, rho * 0.99);

%%
% Step 1 : solve linear system
algo.assign(tmp, KKT \ (rhsKKTconst + rhsKKTvar*(lam-y)));

% Compute over-relaxation
algo.assign(Lx_hat, L*tmp(1:n));

% Step 2 : compute prox operators
algo.assign(q, Lx_hat + l + lam);
prev_y   = y;
for i = 1:nProx
 ind = proxInd{i};
 y(ind) = prox(i).func(q(ind));
end

% Update the lagrange multiplier
%lam = q - y;

% Test: apply more dual steps than primal steps, say 3 dual steps and 1 primal steps
for j= 1:1
 lam = Lx_hat + l + lam - y;
end

% Convergence checks
s = -rho*sum_L*(y - prev_y);
r = Lx_hat + l - y;

rDual   = norm(s);
rPrimal = norm(r);

if rDual < settings.dualTol && rPrimal < settings.primalTol
 cprintf([0.1 0.5 0.1], '\n\n>>>>> Stopping on optimality after %i steps <<<<<\n\n', iii);
 break
end
