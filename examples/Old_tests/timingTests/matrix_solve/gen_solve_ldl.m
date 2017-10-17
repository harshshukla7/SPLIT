function gen_solve_ldl(A,Astr,funcName)

% Pre-factor the matrix A
A = sparse(A);
[L,D,p,S] = ldl(A,0.01,'vector')
s = diag(S);
s = s(p);

% y  = L \ (s .* b(p));
lower_triangular_solve(dat, L, Lstr, s, p, funcname)

% z  = inv(D) * y;
gen.ge

q  = L' \ z;
x(p,1)  = s .* q;


clear

n = 10;
A = [];
while rank(full(A)) < n
  A = sprandn(n,n,0.2);
  A = triu(A);
  A = A + A';
end

[L,D,p,S] = ldl(A,0.01,'vector')
s = diag(S);
s = s(p);

b = randn(size(A,1),1);

y  = L \ (s .* b(p));
z  = D \ y;
q  = L' \ z;
x(p,1)  = s .* q;

sol = inv(A)*b;

[full(x) sol]

% L  * y = P'*S*b
% D  * z = y
% L' * q = z
% x = S * P * q


% 
% 
% x(p,:) = Lm'\(Dm\(Lm\(bM(pm,:))));
% 
% ===
% P'*S*A*S*P = L*D*L'
% ===
% 
% A*x = b
% 
% inv(S)*P*L*D*L'*P'*inv(S) * x = b
% 
% L*D*L'*P'*inv(S) * x = P'*S*b
% 
% ------
% L  * y = P'*S*b
% D  * z = y
% L' * q = z
% x = S * P * q
% 
% 
% 
% M = A;
% [Lm, Dm, pm] = ldl(M, 'vector');
% fprintf(1, 'The error norm ||M(pm,pm) - Lm*Dm*Lm''|| is %g\n', ...
%   norm(full(M(pm,pm) - Lm*Dm*Lm')));
% 
% bM = randn(size(M,1),1);
% sol = inv(A)*bM;
% 
% % Solve a system with this kind of factorization.
% clear x;
% x(pm,:) = Lm'\(Dm\(Lm\(bM(pm,:))));
% fprintf('The absolute error norm ||x - sol|| is %g\n', ...
%   norm(x - sol));
