%% Generate random problem

clear all
close all


n = 2;
m = 1;
N = 50;

[A,B,C,D] = ssdata(drss(n,n,m));

x0 = 10*ones(n,1);


%% Define weight matrix

Q = eye(n);
R = eye(m);
S = zeros(m,n);
[~,P,~]=dlqr(A,B,Q,R);

%% Try with quadprog problem
%%%% Formation of the problem is in following form 
%%%%% Aeq [x;u] = beq
%%%%%% Aeq has the following form 
%%%%%% -I(n)   ..     ..     ..  B
%%%%%%  A     -I(n)   ..     ..  ..  B
%%%%%%  0       A    -I(n)   ..  ..  .. B

%%%%%% beq = zeros except first n elements which are -Ax


H = blkdiag(kron(eye(N),Q),kron(eye(N),R));
H((N-1)*n+1:N*n,(N-1)*n+1:N*n) = P;

%%% next stack for the equality constraints

Aeq = sparse((N)*(n),(N)*(n+m));
beq = sparse((N)*(n),1);

beq(1:n) = -A*x0;
Aeq(1:n,1:n) = -speye(n);

for i=1:N-1
    
    Aeq((i*n+1:(i+1)*n),(i-1)*n+1:i*n) = A;
    Aeq(i*n+1:(i+1)*n,i*n+1:(i+1)*n) = -speye(n);
    Aeq((i-1)*n+1:(i)*n,(n)*N+(i-1)*m+1:n*N+i*m) = B;
    
end

Aeq((N-1)*n+1:(N)*n,(n)*N+(N-1)*m+1:n*N+N*m) = B;

sol = quadprog(H,[],[],[],Aeq,beq);


%% solve using riccati recursion


[K_cont,P] = riccati_recurison_quad(A,B,P,Q,R,S,N);
        
x = x0;

for i=0:N-1
    
    u = K_cont{1,i+1}*x;
    x = A*x + B*u;

end

%% compare the results



%% plot the sparsity of the matrix




%% SPLIT

% x = splitvar(n,N);
% u = splitvar(m,N-1);
%
% x0 = parameter(n,1);
%
% %r = parameter(n,N);
% x(:,1) == x0;
% x(:,2:end) == A*x(:,1:end-1) + B*u;
%
% obj = 0;
% for i=1:N-1
%     obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i);
% end
% obj = obj + x(:,N)'*Q*x(:,N);
%
%
% minimize(obj);
%
% prob = splitProb.genProblem;
%
% par_ex = 10*ones(n,1);
% x0.set(par_ex);
%
% settings.adaptive = 'no';
% settings.precond = 'no';
% sd = coder_admm(prob,settings);
% sol    = AdPrADMM(prob);
% sol_x  = sol.x;
%
% sd.add_var('par_ex', par_ex);
% sd.add_var('sol_x',  sol_x);


