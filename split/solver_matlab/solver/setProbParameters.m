function prob = setProbParameters(paramProb)
%
% prob = setParameters(paramProb)
%
% Compute a static version of the problem by setting all the parameters to
% their current values
%

param = splitProb.getVarValue(paramProb.param);

% Convert the equality constraints
pDat = paramProb.dat;
dat.A = pDat.A;
dat.b = pDat.b + pDat.pB*param;

% Convert the objective
%  0.5*x'*Q*x + (pF*y + f)'*x + (0.5*y'*pQ*y + pf'*y + c)
dat.Q = pDat.Q;
dat.f = param'*pDat.pF + pDat.f;
dat.c = pDat.c + pDat.pf*param + 0.5*param'*pDat.pQ*param;

% Convert the prox functions
prox = paramProb.prox;
for i = 1:length(prox)
  prox(i).l = prox(i).l + prox(i).pL*param;
end
prox = rmfield(prox, 'pL');

prob.prox   = prox;
prob.dat    = dat;
prob.vars   = paramProb.vars;
