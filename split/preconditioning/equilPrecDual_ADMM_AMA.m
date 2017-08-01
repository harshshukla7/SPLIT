function [dualPrec, rho] = equilPrecDual_ADMM_AMA(prob, flagMthd, type)
% Returns dual preconditioner E and optimal stepsize rho

H = full(prob.dat.Q);
A = [];
for iii = 1:length(prob.prox)
    A = [A;prob.prox(iii).L];
end
%% 
switch type
    case 'sdp'
        if (min(eig(H))-1e-6 >= 0) || (isempty(prob.dat.A))
            cvx_begin sdp
            cvx_solver mosek
            variable invP(size(A,1), size(A,1)) diagonal
            variable t
            minimize ( t )
            subject to
            t*H >= invP;
            H <= invP;
            cvx_end
            P = inv(full(invP));
            E = P^(1/2);
            dualPrec = E;
            if strcmp(flagMthd,'fama')
                E = E./sqrt((max(eig(E*A*(H\A')*E))));
            end
            dualPrec = E; 
        end
        if (min(eig(H))-1e-6 < 0) && (~isempty(prob.dat.A))
            K = inv([H prob.dat.A'; prob.dat.A zeros(size(prob.dat.A,1),size(prob.dat.A,1))]);
            K11 = K(1:size(H,1),1:size(H,1));
        end
    case 'equilibration'
        % Equilibration Knight11'
        Q = A*(H\A');
        D = diag(ones(size(A,1),1));  
        E = diag(ones(size(A,1),1)); 
        r = zeros(size(A,1),1);
        c = zeros(size(A,1),1);
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

        if strcmp(flagMthd,'fama')
            E = E./sqrt((max(eig(E*A*(H\A')*E))));
        end
        dualPrec = E;
end
    
%%
% Parameter selection for ADMM and FAMA
if strcmp(flagMthd,'admm')
        lamb_FHF = svd(dualPrec*A*(H\(A'*dualPrec')));
        eig_max  = max(lamb_FHF);
        temp     = find(lamb_FHF>1e-8);
        eig_min  = min(lamb_FHF(temp));
        rho      = 1/sqrt(eig_max*eig_min);
elseif strcmp(flagMthd,'fama') 
    rho = min(eig(H)) / max(eig(A'*A));
end
    
end