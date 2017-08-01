function [rho, alpha, P2, status] = optimal_selections_linear_ADMM(F,H,f)
%>>>> Works for strongly convex QPs only <<<<%

% Optimal stepsize and relaxation parameters
lamb_FHF = svd(F*(H\F'));
eig_max  = max(lamb_FHF);
temp     = find(lamb_FHF>1e-6);
eig_min  = min(lamb_FHF(temp));
rho      = 1/sqrt(eig_max*eig_min);
alpha    = 2;

% Optimal diagonal preconditioning
NF = null(F');
cvx_begin sdp
    cvx_solver mosek
    variable P2(size(F,1),size(F,1)) diagonal
    variable s(1,1)
    minimize ( -s )
    subject to
    diag(P2) >= eps;
    P2 - s*F*F' >= 0;
    NF'*(F*F'-P2)*NF >= 0;
cvx_end
P2 = full(P2);

if norm(P2,'fro') <= 1e-5
    disp('Sd ill-conditioned: ignore preconditioning')
    status = 'bad';
else
    F = P2*F;   f=P2*f;
    lamb_FHF = svd(F*(H\F'));
    eig_max  = max(lamb_FHF);
    temp     = find(lamb_FHF>1e-8);
    eig_min  = min(lamb_FHF(temp));
    rho      = 1/sqrt(eig_max*eig_min);
end
    