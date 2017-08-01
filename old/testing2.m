clear all;

N = 10;
n = 3;
m = 2;

sys = drss(n,1,m);
A = sys.a;
B = sys.b;
C = sys.c;
D = sys.d;
Fx = randn(10,n);
fx = abs(randn(10,1));
Fu = randn(5,m);
fu = abs(randn(5,1));
x0 = randn(n,1);

Q = randn(n,n);
Q = Q'*Q;
R = eye(m);



cvx_begin
cvx_solver sdpt3
variables x(n,N) u(m,N-1)
obj = 0;
for t = 1:N-1
    obj = obj + x(:,t)'*Q*x(:,t) + u(:,t)'*R*u(:,t);
end
obj = obj + x(:,N)'*Q*x(:,N);
minimize ( obj )
subject to
for t = 2:N-1
    Fx*(x(:,t)) <= fx;
    Fu*(u(:,t)) <= fu;
end
x(:,2:N) == A*x(:,1:N-1) + B*u;
x(:,1) == x0;
Fx*(x(:,N)) <= fx;
cvx_end

splitProb.clearProblem;
X = splitvar(n,N);
U = splitvar(m,N-1);
obj = 0;
for t = 1:N-1
    obj = obj + X(:,t)'*Q*X(:,t) + U(:,t)'*R*U(:,t);
end
obj = obj + X(:,N)'*Q*X(:,N);
for t = 2:N-1
    Fx*(X(:,t)) <= fx;
    Fu*(U(:,t)) <= fu;
end
X(:,2:N) == A*X(:,1:N-1) + B*U;
X(:,1) == x0;
Fx*(X(:,N)) <= fx;

minimize(obj);

prob = splitProb.genProblem;
settings.maxItr   = 3e3;
settings.dualTol = 1e-2;
settings.primalTol = 1e-2;
        
% w/o rstart    
[sol, stats, stats_plot] = PrFAMA(prob, settings);
X = reshape(sol.x(1:n*N),n,N);
U = reshape(sol.x(n*N+1:end),m,N-1);

% w/ restart
settings.restart = 'yes';
[sol, stats, stats_plot] = PrFAMA(prob, settings);
X = reshape(sol.x(1:n*N),n,N);
U = reshape(sol.x(n*N+1:end),m,N-1);

