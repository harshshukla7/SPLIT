function  genProx(cdr, q, rho, y)
%
% Call the proximal function
%
%  y = prox(q, rho)
%
% y, q and rho are strings describing these vectors
%
% Allocates memory for the vector y


cdr.print('%s = proxFunc(%s, %s, %s);\n', y, q, rho, y)
cdr.addVar(y, length(cdr.l), 1);

% Queue the creation of the proximal function
cdr.addFunction('function y = proxFunc(q, rho, y)',...
  cdr.genProxFunction('q','y','rho'),...
  sprintf('Compute the proximal funciton'))
