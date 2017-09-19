%% PDOS diagonal preconditioning

function [D,E] = PDOS_diag_prec(T,mi,M)
% function [] = PDOS_diag_prec(m,n,mi,M)
% ----------------------------------------------------------
% ----------------------------------------------------------
% Generates diagonal preconditioners D, E for the problem
% minimize 1/2*z'*Q*z + sum_M li(Ti*z+ti)
% subject to A*z = b,
% such that the matrix D*[T1; T2; ... TM; A]*E has a smaller spectral
% radius.
% D = diag(Di), where Di is of size mi x n and sum_M mi = m.
% ----------------------------------------------------------
% ----------------------------------------------------------

n = size(T,2);
Ttilde = zeros(M,n);
mi = [0; mi];
for iii = 2:M+1
    Ttilde(iii-1,:) = (1/mi(iii)) * sum(T(sum(mi(1:iii-1))+1:sum(mi(1:iii)),:).^2,1);  
end

p = zeros(M,1);
for iii = 1:M
    for jjj = 1:n
        p(iii) = p(iii) + Ttilde(iii,jjj);
    end
    p(iii) = 1 / sqrt(p(iii));
end

delta = zeros(n,1);
for jjj = 1:n
    for iii = 1:M
        delta(jjj) = delta(jjj) + Ttilde(iii,jjj) * p(iii)^2;
    end
    delta(jjj) = 1 / sqrt(delta(jjj));
end

D = [];
for iii = 1:M
    D = blkdiag(D, p(iii)*eye(mi(iii+1)));
end
E = diag(delta);
    
    
