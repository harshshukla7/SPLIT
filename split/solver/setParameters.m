function prob = setProbParameters(paramProb)
%
% prob = setParameters(paramProb)
%
% Compute a static version of the problem by setting all the parameters to
% their current values
%

dat    = paramProb.dat;


prob.prox   = paramProb.prox;
prob.varInd = paramProb.varInd;
