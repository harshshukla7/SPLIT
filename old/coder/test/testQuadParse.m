splitProb.clearProblem

% Testing the parsing of the problem
%
%  h'*x + x'*Q*x + h'*x
%

Q = eye(2);
h = [1;1];

x   = splitvar(2,1);
par = parameter(3,1);
F   = 3*ones(3,2);

obj = h'*x + x'*Q*x + par'*F*x;

minimize(obj);

prob = splitProb.genProblem
