function [sumMu,Tau,Sigma] = equilPrecLeftRightPDA(prob, alpha, precDual, precPrimal)

%%
% Function description

%%

K = [];
for iii = 1:length(prob.prox)
    K = [K;prob.prox(iii).L];
end
[m,n] = size(K);
if any(prob.dat.Q)
    Lip = max(eig(prob.dat.Q));
else
    Lip = 0;
end

%% only primal (dual is given)
% given dual corresponds to precDual=E=Sigma^(1/2)
if ~isempty(precDual)
    alpha = 0;
    % dual preconditioner is provided; find primal that ensures bounded unit norm
    Kd = precDual*K; % Kd = Sigma^(1/2)*K
    sigma = ones(m,1);
    sum2 = zeros(n,1);
    for jjj = 1:n
        for iii = 1:m
            sum2(jjj) = sum2(jjj) + abs(Kd(iii,jjj))^(2-alpha);
        end
        tau(jjj) = 1/sum2(jjj);
    end
    Sigma = diag(sigma);    Tau = diag(tau);
    Mat = Sigma^(1/2)*Kd*Tau^(1/2); % Sigma is identity since dual preconditioner already included in K
    Sigma = precDual^2; % recover dual as provided
    Mat = Sigma^(1/2)*K*Tau^(1/2);
    
    sumMu = 0;  ni_old = 0;
    for iii = 1:length(prob.prox)
        ni = size(prob.prox(iii).L,1);
        Pi = precDual(ni_old+1:ni_old+ni,ni_old+1:ni_old+ni);
        sumMu = sumMu + max(abs(eig((Pi*prob.prox(iii).L*Tau^(1/2))'*(Pi*prob.prox(iii).L*Tau^(1/2)))));
        ni_old = ni;
    end
    while ( sqrt(sumMu) >= 1-(Lip/2)*sqrt(max(abs(eig(Tau'*Tau)))) )
        Tau = Tau.*0.99;
        sumMu = 0;  ni_old = 0;
        for iii = 1:length(prob.prox)
            ni = size(prob.prox(iii).L,1);
            Pi = precDual(ni_old+1:ni_old+ni,ni_old+1:ni_old+ni);
            sumMu = sumMu + max(abs(eig((Pi*prob.prox(iii).L*Tau^(1/2))'*(Pi*prob.prox(iii).L*Tau^(1/2)))));
            ni_old = ni;
        end
    end
    Mat = Sigma^(1/2)*K*Tau^(1/2);
    
    sprintf('condition number is %4.4f', max(abs(eig(Mat'*Mat))))
end

%% only dual (primal is given)
% given dual corresponds to precDual=E=Sigma^(1/2)
if ~isempty(precPrimal)
    alpha = 2;
    % dual preconditioner is provided; find primal that ensures bounded unit norm
    Kp = K*precPrimal; % Kp = K*Tau^(1/2)
    tau = ones(n,1);
    sum1 = zeros(m,1);
    for iii = 1:m
        for jjj = 1:n
            sum1(iii) = sum1(iii) + abs(Kp(iii,jjj))^alpha;
        end
        sigma(iii) = 1/sum1(iii);
    end 
    Sigma = diag(sigma);    Tau = diag(tau);
    Mat = Sigma^(1/2)*Kp*Tau^(1/2); % Tau is identity since dual preconditioner already included in K
    Sigma = Sigma./(max(eig(Mat'*Mat))); % primal scaled the norm
    Tau = precPrimal^2; % recover dual as provided
    precDual = Sigma^(1/2);
    Mat = Sigma^(1/2)*K*Tau^(1/2);
    
    sumMu = 0;  ni_old = 0;
    for iii = 1:length(prob.prox)
        ni = size(prob.prox(iii).L,1);
        Pi = precDual(ni_old+1:ni_old+ni,ni_old+1:ni_old+ni);
        sumMu = sumMu + max(abs(eig((Pi*prob.prox(iii).L*Tau^(1/2))'*(Pi*prob.prox(iii).L*Tau^(1/2)))));
        ni_old = ni;
    end
    while ( sqrt(sumMu) >= 1-(Lip/2)*sqrt(max(abs(eig(Tau'*Tau)))) )
        Tau = Tau.*0.99;
        sumMu = 0;  ni_old = 0;
        for iii = 1:length(prob.prox)
            ni = size(prob.prox(iii).L,1);
            Pi = precDual(ni_old+1:ni_old+ni,ni_old+1:ni_old+ni);
            sumMu = sumMu + max(abs(eig((Pi*prob.prox(iii).L*Tau^(1/2))'*(Pi*prob.prox(iii).L*Tau^(1/2)))));
            ni_old = ni;
        end
    end
    Mat = Sigma^(1/2)*K*Tau^(1/2);
    
    sprintf('condition number is %4.4f', sqrt(max(abs(eig(Mat'*Mat))))^2)
end

%% perform full equilibration
if isempty(precDual) && isempty(precPrimal)
    tau = ones(n,1);    sigma = ones(m,1);
    sum1 = zeros(m,1);
    for iii = 1:m
        for jjj = 1:n
            sum1(iii) = sum1(iii) + abs(K(iii,jjj))^alpha;
        end
        sigma(iii) = 1/sum1(iii);
    end    
    sum2 = zeros(n,1);
    for jjj = 1:n
        for iii = 1:m
            sum2(jjj) = sum2(jjj) + abs(K(iii,jjj))^(2-alpha);
        end
        tau(jjj) = 1/sum2(jjj);
    end
    Sigma = diag(sigma);    Tau = diag(tau);
    Mat = Sigma^(1/2)*K*Tau^(1/2);
    precDual = Sigma^(1/2);

    sumMu = 0;  ni_old = 0;
    for iii = 1:length(prob.prox)
        ni = size(prob.prox(iii).L,1);
        Pi = precDual(ni_old+1:ni_old+ni,ni_old+1:ni_old+ni);
        sumMu = sumMu + max(abs(eig((Pi*prob.prox(iii).L*Tau^(1/2))'*(Pi*prob.prox(iii).L*Tau^(1/2)))));
        ni_old = ni;
    end
    while ( sqrt(sumMu) >= 1-(Lip/2)*sqrt(max(abs(eig(Tau'*Tau)))) )
        Tau = Tau.*0.99;
        sumMu = 0;  ni_old = 0;
        for iii = 1:length(prob.prox)
            ni = size(prob.prox(iii).L,1);
            Pi = precDual(ni_old+1:ni_old+ni,ni_old+1:ni_old+ni);
            sumMu = sumMu + max(abs(eig((Pi*prob.prox(iii).L*Tau^(1/2))'*(Pi*prob.prox(iii).L*Tau^(1/2)))));
            ni_old = ni;
        end
    end
    Mat = Sigma^(1/2)*K*Tau^(1/2);
    
    sprintf('condition number is %4.4f', max(abs(eig(Mat'*Mat))))
end 
    

end
