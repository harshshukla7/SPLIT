clear

splitProb.clearProblem;

x = splitvar(3,1);
y = parameter(2,1);

Q = randn(3); Q=Q*Q';
F = randn(3,2);

A = randn(20,3);
A*x <= 10*ones(20,1);

B = randn(1,3);
B*x == 1;

% norm(randn(3,3)*x + randn(3,1),inf) <= sum(x);
minimize(y'*F'*x);% + x'*Q*x);% + 3.4*norm(x))

prob = splitProb.genProblem;

xx = randn(2,1);



%%

cdr = splitCoder(prob, 'admm', 'funcADMM')

cdr.rho = 1;
cdr.KKTSolveMode = 'Invert KKT';
% cdr.KKTSolveMode = 'LU';
cdr.MatVec_SparseLimit = .1;
cdr.itr_conv = 1;

solver = cdr.genMatlab;


[x, y, lambda, itr, tm] = solver(xx);

yal_x = sdpvar(3,1);
yal_y = xx;

con = set(A*yal_x <= 10*ones(20,1));
% obj = yal_y'*F'*yal_x + yal_x'*Q*yal_x + 3.4*norm(yal_x);
obj = yal_y'*F'*yal_x;

solvesdp(con, obj)


cdr.genMex;

[x, y, lambda, itr, tm] = funcADMM_mex(xx);


% pause

%%

cdr = splitCoder2(prob, 'template_cpii.m', 'funcCPII')

% cdr.MatVec_SparseLimit = 0.1;
cdr.itr_conv = 10;
cdr.maxItr = 1e6;

solver = cdr.genMatlab;

[x,y,itr,tm] = solver(xx);

cdr.genMex;

x2 = x;

% pause


%%

% 
% N = 1e4;
% fprintf('Solving %i random problems...\n', N);
% PAR = randn(2,N);
% tic
% for i = 1:N
%   [x,rDual,rPrimal,itr(i),tm(i)] = funcADMM_mex(PAR(:,i));
% end
% t = toc;
% 
% fprintf('Time per iteration = %.2e\n', t / sum(itr));
% fprintf('Time per solve     = %.2e\n', t / N);
