
clear all
%% parameters to play with
n = 3;
m = 3;
N = 6;


rho = 1;
[A,B,C,D] = ssdata(drss(n,n,m));

x = splitvar(n,N);
u = splitvar(m,N-1);
x0 = parameter(n,1);

constraint = 'box';

switch constraint
    
    case 'box'
    -10 <= x(:) <= 10;
    -1 <= u(:) <= 1;

    case 'l1_ball'
        
        %%%%% 2. l1-ball
    F = 9e-3*randn(n*N,n);
    f = 1e-6*randn(n,1);
    d = 100;
    norm(F'*x(:)+f,1) <= d;
    
    case 'l2_ball'
    %%%%% 2. l2-ball
    F = 50*randn(n*N,n);
    f = 1e-7*randn(n,1);
    d = 1;
    norm(F'*x(:)+f,2) <= d;
    
    case 'linf_ball'
    %%%%% 3. linf ball
    F = 15*randn(n*N,n);
    f = 1e-5*randn(n,1);
    d = 1;
    norm(F'*x(:)+f,inf) <= d;
    
    case 'lorentz_cone'
    %%%%% 4. Lorentz Cone

    F = 50*randn(n*N,n);
    f = 1e-7*randn(n,1);
    d = 1;
    q = 1e-3*randn(n*N,1);

    norm(F'*x(:)+f,2) <= q'*x(:) + d;

    case 'polytope'
    %%%%% 5. Polytope
    K = 1e-1*randn(3*(n*N+m*(N-1)),n*N+m*(N-1));
    k = 20*abs(randn(3*(n*N+m*(N-1)),1));
    K*[x(:); u(:)] <= k;
    
end 


x(:,1) == x0;

x(:,2:end) == A*x(:,1:end-1) + B*u;

pattern = 'diagonal'; %%%% there are 3 cases

switch pattern
    
    case 'diagonal'
    %% Q diagonal

    Q = eye(n);

    case 'banded'
    %%  Q banded

    band = 3;
    Q = createbandOPT(n,band);

    
    case 'full'
    %% Q full

    Q = ones(n,n);

end
R = eye(m);
objX = x .* (Q*x);
objU = u .* (R*u);
obj = sum(objX(:)) + sum(objU(:));

minimize(obj);

prob = splitProb.genProblem;
dat = prob.coderData;

Q   = dat.Q;
L   = dat.L;
A   = dat.A;
K = [Q + rho*L'*L, A'; A, zeros(size(A,1))];
K1  = [Q , A'; A, zeros(size(A,1))];

figure()
spy(K)

figure()
spy(K1)


%% ldl decomposition

[L,D,T]=ldl(K);

figure()
spy(L)