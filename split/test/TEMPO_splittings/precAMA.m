function [dualPrec] = precAMA(dat)
% Returns dual preconditioner E

% import data
nx = dat.dim.nx;
nu = dat.dim.nu;
N = dat.time.N;
Q = dat.cost.Q;
R = dat.cost.R;
P = dat.cost.P;
L = dat.L;

H = 2*blkdiag(kron(eye(N),Q),P,kron(eye(N),R));
%% 
% Equilibration Knight11'
Q = L*(H\L');
D = diag(ones(size(L,1),1));  
E = diag(ones(size(L,1),1)); 
r = zeros(size(L,1),1);
c = zeros(size(L,1),1);
for kkk = 1:1e1
    for i = 1:size(Q,1)
        r(i) = sqrt(norm(Q(i,:),Inf));
    end
    R = diag(r);
    for j = 1:size(Q,2)
        c(j) = sqrt(norm(Q(:,j),Inf));
    end
    C = diag(c);
    Q = R\(Q*inv(C));
    D = D*inv(R);
    E = E*inv(C);
end

E = E./sqrt((max(eig(E*L*(H\L')*E))));
dualPrec = E;
end