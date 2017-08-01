function [Tau, Sigma] = CP_diag_prec(K)

 alpha = 1;
 sumsigma = zeros(size(K,1),1);
 sumtau = zeros(size(K,2),1);
 for j = 1:size(K,2)
     for i = 1:size(K,1)
         sumtau(j) = sumtau(j) + abs(K(i,j))^(2-alpha);
         sumsigma(i) = sumsigma(i) + abs(K(i,j))^alpha;
     end
 end

 Tau = diag(1./sumtau);
 Sigma = diag(1./sumsigma);

