function [dualPrec, rho] = optPrecDual_ADMM_AMA(prob, flagMthd)
%% 
% SDP preconditioning
C = full(prob.dat.A);
H = full(prob.dat.Q);
A = full(prob.prox.L);
if (min(eig(A*(H\A')))>eps) % positive Hessian
   cvx_begin sdp
   cvx_solver sdpt3
   variable M(size(A,1),size(A,1)) diagonal
   variable s(1,1)
   minimize ( s )
   subject to
   s*(A*(H\A')) >= M;
   (A*(H\A')) <= M;
   cvx_end  
   M = full(M);
   if norm(M,'fro') <= 1e-5
        disp('Sd ill-conditioned: ignore preconditioning')
        status = 'bad';
   end
   dualPrec = (inv(full(M)))^(1/2);
elseif (min(max(0,eig(A*(H\A'))))==0) && ~isempty(C) % zero-eig Hessian 
    K = inv([H C'; C zeros(size(C,1),size(C,1))]);
    [V,D] = eig(full(K(1:size(H,1),1:size(H,1))));
    K11 = V*max(D,0)*V';
    [V,D] = eig(A*K11*A');
    R = D^(1/2)*V';
    cvx_begin sdp
    cvx_solver mosek
    variable M(size(R,2),size(R,2)) diagonal
    variable s(1,1)
    minimize ( -s )
    subject to
    R*M*R' <= eye(size(R,1));
    R*M*R' >= s*eye(size(R,1));
    cvx_end
    M = full(M);
    if norm(M,'fro') <= 1e-5
        disp('Sd ill-conditioned: ignore preconditioning')
        status = 'bad';
    end
    dualPrec = full(M)^(1/2);
elseif (min(max(0,eig(A*(H\A'))))==0) && isempty(C)
    [rho, alpha, dualPrec, status] = optimal_selections_linear_ADMM(A,H);
end

if strcmp(status,'bad')
    disp('Ill-conditioned preconditioner after solving the SDP - try with equilibration')
end

%%
% Parameter selection for ADMM and FAMA
if strcmp(flagMthd,'admm')
    if ~strcmp(status, 'bad')
        lamb_FHF = svd(dualPrec*A*(H\(A'*dualPrec')));
        eig_max  = max(lamb_FHF);
        temp     = find(lamb_FHF>1e-8);
        eig_min  = min(lamb_FHF(temp));
        rho      = 1/sqrt(eig_max*eig_min);
    else
        lamb_FHF = svd(A*(H\A'));
        eig_max  = max(lamb_FHF);
        temp     = find(lamb_FHF>1e-8);
        eig_min  = min(lamb_FHF(temp));
        rho      = 1/sqrt(eig_max*eig_min);
    end
elseif strcmp(flagMthd,'ama') && strcmp(status, 'bad')
    rho = min(eig(H)) / max(eig(A'*A));
end
end
   