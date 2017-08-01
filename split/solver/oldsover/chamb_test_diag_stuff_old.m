% number of elements
n = 10;

% parameter
tauprev = 0.985;

% generate a random vector bro:
v = rand(n,1)*3.25 - 0.02*rand(n,1)*6.5;

% diagonal matrix plus and identity stuff
Q = diag(v);
I = eye(n);

% expression to be inverted
expression = Q + 1/tauprev*I;

% compute inverse by the general approach
invExpr = inv(expression);

% compute inverse by the stuff
hab = v + diag(I)*1/tauprev;
hab = 1./hab;



