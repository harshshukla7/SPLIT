
function sd = coder_PDA(prob,settings)
% Generate custom c-code to solve PDA
%% settings.E is such that yd = Ey where y is unscaled variables
%  settings.Dinv is such that and Dinv = inv(D) and Zp = Dz where z is unscaled variables


%% Intial Checks

warning('Dimensions of matrices and vectors provided by user are not checked.');

if(~isfield(settings, 'precond'))
    settings.precond = 'no';
end

if(~isfield(settings, 'adaptive'))
    settings.adaptive = 'no';
end

if(strcmp(settings.precond,'yes') && strcmp(settings.precond,'no'))
    error('precond value must be either "yes" or "no"');
end

if(strcmp(settings.adaptive,'yes') && strcmp(settings.adaptive,'no'))
    error('adaptive value must be either "yes" or "no"');
end


if(strcmp(settings.precond, 'yes') && ~isfield(settings, 'scale_rDual'))
    error('User must provide scaling Dual scaling for preconditioning')
end

if(strcmp(settings.precond, 'yes') && ~isfield(settings, 'scale_rPrimal'))
    error('User must provide scaling matrix Primal scaling for preconditioning')
end

if(isfield(settings, 'scale_rDual') && ~isdiag(settings.scale_rDual))
    error('Scaling matrix for Dual must be diagonal')
end

if(isfield(settings, 'scale_rPrimal') && ~isdiag(settings.scale_rPrimal))
    error('Scaling matrix D_scale must be diagonal')
end

if(strcmp(settings.precond, 'yes') && ~isfield(settings, 'Mu'))
    error('User must provide Mu for preconditioning')
end


%%%%% next set alpha gamma

if ~isfield(settings, 'relaxation'), settings.relaxation = 1; end

if(isfield(settings, 'Mat_Vec'))
    Mat_Vec = settings.Mat_Vec;
end

if(~isfield(settings, 'Mat_Vec'))
    Mat_Vec = 'auto';
end

if(isfield(settings, 'Lin_Solve'))
    Lin_Solve = settings.Lin_Solve;
end

if(~isfield(settings, 'Lin_Solve'))
    Lin_Solve = 'auto';
end



%% Extrace the data from SpliProb

dat  = prob.coderData;
prox  = prob.prox;
nProx = length(prox);

sd   = splitData;
dat.A = sparse(dat.A);
R     = chol(dat.A*dat.A');


if(isfield(settings, 'primal_vars_x'))
    %% ifdef warm start x
    sd.hl('#define warm_start_x \n');
    a_name=sprintf('warm_x');
    sd.add_var(a_name,settings.primal_vars_x, 'type', 'real');
    
    
end

% if(isfield(settings, 'primal_vars_y'))
%      %% ifdef warm start y
%      sd.hl('#define warm_start_y \n');
%      a_name=sprintf('warm_y');
%     sd.add_var(a_name,settings.primal_vars_y, 'type', 'real');
% end

if(isfield(settings, 'dual_vars_lam'))
    %% ifdef warm start lam
    sd.hl('#define warm_start_lam \n');
    a_name=sprintf('warm_lambda');
    
    if(strcmp(settings.precond,'yes'))
        sd.add_var(a_name,settings.E*settings.dual_vars_lam, 'type', 'real');
    else
        sd.add_var(a_name,settings.dual_vars_lam, 'type', 'real');
    end
    
end


if(strcmp(settings.adaptive,'yes'))
    sd.hl('#define adaptive\n');
end


if(isfield(settings, 'adaptive_everystep'))
    if(strcmp(settings.adaptive_everystep,'yes'))
        sd.hl('#define adaptive_everystep\n');
    end
end


if(strcmp(settings.precond,'yes'))
    sd.hl('#define precond\n');
end

sd.define('nParam',  size(dat.pL,2), 'int');
sd.define('nPrimal', size(dat.A,2),  'int');
sd.define('nDual',   size(dat.L,1),  'int');
sd.define('nEqCon',  size(dat.A,1),  'int');


%% Data
n = size(dat.Q, 1);
m = size(dat.A, 1);

L = []; l = [];
nProxVars = zeros(nProx,1);
kk = 1;
sumMu = 0;
for i = 1:nProx
    sumMu        = sumMu + max(svds(prox(i).L))^2;
    l            = [l;prox(i).l];
    L            = [L;prox(i).L];
    proxInd{i}   = kk + [0:size(prox(i).L,1)-1];
    kk           = kk + size(prox(i).L,1);
    nProxVars(i) = size(prox(i).L,1);
end
proxWeight = [prob.prox.weight];



%% Findout tau,rho,thetha,delta, kappa for not precondtion

if strcmp(settings.precond,'no') % No prec, No adapt
    
    settings.Mu = max(svds(L))^2; % less conservative bound
    Mu = settings.Mu;
    if (~isempty(prob.dat.Q))&&(nnz(prob.dat.Q)) % existence of smooth term changes stepsize
        settings.tau = 1/max(eig(dat.Q));
        tau = settings.tau;
        settings.rho = (1-tau*max(eig(dat.Q))/2)^2 / (tau*Mu);
        rho = settings.rho;
        % over-relaxation
        settings.kappa = (1/max(eig(dat.Q)))*(1/tau-rho*Mu);
        kappa = settings.kappa;
        settings.delta = max(4*kappa/(2*kappa+1),min(3/2,1/2+kappa));
        delta = settings.delta;
        settings.theta = delta;
        theta = delta;
    else
        settings.tau = 1/sqrt(Mu);
        tau = settings.tau;
        settings.rho = 1/(tau*Mu);
        rho = settings.rho;
        % over-relaxation
        settings.theta = 1.8;
    end
    if (strcmp(settings.adaptive,'yes')) % No prec, Yes adapt
        if any(dat.Q)
            settings.tau = 0.05/max(eig(dat.Q));
            tau = settings.tau;
            settings.rho = (1-tau*max(eig(dat.Q))/2)^2 / (tau*sumMu);
            rho = settings.rho;
        else
            settings.tau = 1/sqrt(sumMu);
            tau = settings.tau;
            settings.rho = tau;
            rho = settings.rho;
        end
        settings.theta = 1;
        theta = settings.theta;
        % adaptive stepsize constants
        settings.alpha = 0.5;
        alpha = settings.alpha;
        settings.Delta = 1.5;
        Delta = settings.Delta;
        settings.eta = 0.95;
        eta = settings.eta;
        settings.s = 1;
        s = settings.s;
    end
    
    
    %%%% Now write to c
    if(isfield(settings, 'Mu'))
        sd.hl('extern REAL Mu;\n');
        sd.cl(' REAL Mu = %f;\n', settings.Mu);
    end
    
    if(isfield(settings, 'tau'))
        sd.hl('extern REAL tau;\n');
        sd.cl(' REAL tau = %f;\n', settings.tau);
    end
    
    if(isfield(settings, 'rho'))
        sd.hl('extern REAL rho;\n');
        sd.cl(' REAL rho = %f;\n', settings.rho);
    end
    
    if(isfield(settings, 'kappa'))
        sd.hl('extern REAL kappa;\n');
        sd.cl(' REAL kappa = %f;\n', settings.kappa);
    end
    
    if(isfield(settings, 'delta'))
        sd.hl('extern REAL delta;\n');
        sd.cl(' REAL delta = %f;\n', settings.delta);
    end
    
    if(isfield(settings, 'Delta'))
        sd.hl('extern REAL Delta;\n');
        sd.cl(' REAL Delta = %f;\n', settings.Delta);
    end
    
    if(isfield(settings, 'theta'))
        sd.hl('extern REAL theta;\n');
        sd.cl(' REAL theta = %f;\n', settings.theta);
    end
    
    if(isfield(settings, 'alpha'))
        sd.hl('extern REAL alpha;\n');
        sd.cl(' REAL alpha = %f;\n', settings.alpha);
    end
    
    if(isfield(settings, 'eta'))
        sd.hl('extern REAL eta;\n');
        sd.cl(' REAL eta = %f;\n', settings.eta);
    end
    
    if(isfield(settings, 's'))
        sd.hl('extern REAL s_adapt;\n');
        sd.cl(' REAL s_adapt = %f;\n', settings.s);
    end
end



%% if it is preconditoned

if strcmp(settings.precond,'yes')
    
    prec.Dinv = settings.scale_rPrimal; % D^(-1) (preconditioner)
    prec.D = inv(settings.scale_rPrimal); % D
    prec.T = settings.scale_rPrimal*settings.scale_rPrimal'; % T = D^(-1)*D^(-T)
    prec.P  = settings.scale_rDual'*settings.scale_rDual; % P=E'E
    prec.E  = settings.scale_rDual; % E (preconditioner)
    prec.L  = prec.E*L*prec.Dinv;
    
    
    Mu    = settings.Mu;
    theta = 1;
    %% Stepsizes and adaptation constants
    if (strcmp(settings.adaptive,'yes')) % Yes prec, Yes adapt
        if any(dat.Q)
            tau = 1/(max(svds(prec.T))*max(eig(dat.Q)));
            rho = (1-tau*max(svds(prec.T))*max(eig(dat.Q))/2)^2 / (tau*Mu);
            
        else
            tau = 1/sqrt(Mu);
            settings.tau = tau;
            rho = tau;
            settings.rho = rho;
        end
        % adaptive stepsize constants
        alpha = 0.5;
        Delta = 2;
        eta = 0.95;
        s = 1;
        
        %%%% write to c for precond-adaptive case
        
        sd.hl('extern REAL Mu;\n');
        sd.cl(' REAL Mu = %f;\n', Mu);
        
        sd.hl('extern REAL sqrt_Mu;\n');
        sd.cl(' REAL sqrt_Mu = %f;\n', sqrt(Mu));
        
        sd.hl('extern REAL theta;\n');
        sd.cl(' REAL theta = %f;\n', theta);
        
        sd.hl('extern REAL tau;\n');
        sd.cl(' REAL tau = %f;\n', tau);
        
        sd.hl('extern REAL rho;\n');
        sd.cl(' REAL rho = %f;\n', rho);
        
        sd.hl('extern REAL alpha;\n');
        sd.cl(' REAL alpha = %f;\n', alpha);
        
        sd.hl('extern REAL Delta;\n');
        sd.cl(' REAL Delta = %f;\n', Delta);
        
        sd.hl('extern REAL eta;\n');
        sd.cl(' REAL eta = %f;\n', eta);
        
        sd.hl('extern REAL s_adapt;\n');
        sd.cl(' REAL s_adapt = %f;\n', s);
        
    elseif (strcmp(settings.adaptive,'no')) % Yes prec, No adapt
        rho = 1;
        settings.rho = rho;
        tau = 1;
        settings.tau = tau;
        
        sd.hl('extern REAL Mu;\n');
        sd.cl(' REAL Mu = %f;\n', Mu);
        
        sd.hl('extern REAL sqrt_Mu;\n');
        sd.cl(' REAL sqrt_Mu = %f;\n', sqrt(Mu));
        
        sd.hl('extern REAL theta;\n');
        sd.cl(' REAL theta = %f;\n', theta);
        
        sd.hl('extern REAL tau;\n');
        sd.cl(' REAL tau = %f;\n', tau);
        
        sd.hl('extern REAL rho;\n');
        sd.cl(' REAL rho = %f;\n', rho);
        
    end
end



%% custom_compute_parametric

% Compute: l = pL*par + l_const, etc

sd.add_function(coderFunc_parametric('custom_compute_parametric', dat));

%%%% We wont use following variables so set them to zero
sd.define('nn_lp',0,'int');
sd.define('LNZ_ss', 0, 'int');
sd.define('ANZ_ss', 0, 'int');
sd.define('N_ss', 0, 'int');
%     A = 0;
%     sd.add_function(coderFunc_prefactor('custom_compute_prefactor',A,'none'));
%

%% step 1 Linear system solve :


if(strcmp(settings.adaptive,'yes'))
    
    if(strcmp(settings.precond,'yes'))
        
        sd.add_function(...
            coderFunc_times('custom_mult_Ltrans', sparse(dat.L'), 'Astr', 'Ltrans', 'method', Mat_Vec));
        
        
        if (~isdiag(prec.T))
            
            error('prec.T is not diagonal');
        end
        
        a_name=sprintf('prec_T');
        sd.add_var(a_name,diag(prec.T), 'type', 'real');
        
        
        
        TQ = sparse(prec.T*dat.Q);
        sd.add_function(...
            coderFunc_times('custom_mult_TQ', sparse(TQ), 'Astr', 'TQ', 'method', Mat_Vec));
        
    elseif(strcmp(settings.precond,'no'))
        
        sd.add_function(...
            coderFunc_times('custom_mult_Ltrans', sparse(dat.L'), 'Astr', 'Ltrans', 'method', Mat_Vec));
        
        sd.add_function(...
            coderFunc_times('custom_mult_Q', sparse(dat.Q), 'Astr', 'dat_Q', 'method', Mat_Vec));
        
    end
    
    
end

if(strcmp(settings.adaptive,'no'))
    
    if(strcmp(settings.precond,'yes'))
        
        sd.add_function(...
            coderFunc_times('custom_mult_Ltrans', sparse(dat.L'), 'Astr', 'Ltrans', 'method', Mat_Vec));
        
        
        if (~isdiag(prec.T))
            
            error('prec.T is not diagonal');
        end
        
        a_name=sprintf('tau_T');
        sd.add_var(a_name,(settings.tau*diag(prec.T)), 'type', 'real');
        
        
        ITQ =sparse(eye(size(dat.Q))-settings.tau*prec.T*dat.Q);
        sd.add_function(...
            coderFunc_times('custom_mult_ITQ', sparse(ITQ), 'Astr', 'ITQ', 'method', Mat_Vec));
        
    elseif(strcmp(settings.precond,'no'))
        
        sd.add_function(...
            coderFunc_times('custom_mult_Ltrans', sparse(dat.L'), 'Astr', 'Ltrans', 'method', Mat_Vec));
        
        sd.add_function(...
            coderFunc_times('custom_mult_ItauQ', sparse(eye(size(dat.Q))-tau*dat.Q), 'Astr', 'dat_ItauQ', 'method', Mat_Vec));
        
    end
    
    
end


%% Dynamics Update



if ~isempty(dat.A) % dynamics update
    
    sd.hl('#define dynamics_update \n');
    
    chol_solve = coderFunc_chol_solve('custom_chol_solve', (dat.A*dat.A'), 'Astr', 'AA','method', Lin_Solve);
    sd.add_function(chol_solve);
    
    if strcmp(chol_solve.desc,'ldl_ss')
        sd.cl('double *Lx_ss;');
        sd.cl('int *Li_ss;');
        sd.hl('#define suitesparse_linsolve\n');
        
        K = dat.A*dat.A';
        LNZ=nnz(K)-nnz(triu(K));
        ANZ=nnz(triu(K));
        N=size(K,2);
        %%%% add the varibles
        
        sd.define('nn_lp',size(K,2),'int');
        sd.define('LNZ_ss', LNZ, 'int');
        sd.define('ANZ_ss', ANZ, 'int');
        sd.define('N_ss', N, 'int');
        
    end
    
    if strcmp(chol_solve.desc,'ldl_lp')
        sd.cl('__CLPK_integer ipiv[nn_lp];');
        
        sd.hl('#define lapack_linsolve\n');
        K = dat.A*dat.A';
        
        sd.define('nn_lp',size(K,2),'int');
    end
    
    sd.add_function(...
        coderFunc_times('custom_mult_Atrans', sparse(dat.A'), 'Astr', 'dat_Atrans', 'method', Mat_Vec));
    
    %      sd.add_function(...
    %             coderFunc_times('custom_mult_AAA', sparse(dat.A'*(inv(dat.A*dat.A'))), 'Astr', 'AAA', 'method', Mat_Vec));
    %
    
    sd.add_function(...
        coderFunc_times('custom_mult_A', sparse(dat.A), 'Astr', 'dat_A', 'method', Mat_Vec));
    
    
end

%% prox step

sd.hl('extern REAL prox_var; \n');
sd.cl(' REAL prox_var = %f;', 0.0);

sd.add_function(...
    coderFunc_times('custom_mult_L', sparse(dat.L), 'Astr', 'L', 'method', Mat_Vec));

if (strcmp(settings.precond,'yes'))
    
    if (~isdiag(prec.P))
        
        error('prec.P is not diagonal');
    end
    
    a_name=sprintf('prec_P');
    sd.add_var(a_name,diag(prec.P), 'type', 'real');
    
end


% Evaluate prox functions y = prox(workDual)
sd.add_function(coderFunc_prox_conj(dat));


%% residual check

if(strcmp(settings.adaptive,'no'))
    
    if (strcmp(settings.precond,'yes'))
        sd.add_function(...
            coderFunc_times('custom_mult_QTD', sparse(dat.Q-(1/tau)*prec.D), 'Astr', 'QTD', 'method', Mat_Vec));
    end
    
    if (strcmp(settings.precond,'no'))
        sd.add_function(...
            coderFunc_times('custom_mult_QTD', sparse(dat.Q-(1/tau)*eye(size(dat.Q,1),size(dat.Q,2))), 'Astr', 'QTD', 'method', Mat_Vec));
    end
    
end

if (strcmp(settings.precond,'yes'))
    
    
    if (~isdiag(prec.D))
        
        error('prec.D is not diagonal');
    end
    
    a_name=sprintf('prec_D');
    sd.add_var(a_name,diag(prec.D), 'type', 'real');
    
    
    if (~isdiag(prec.Dinv))
        
        error('prec.Dinv is not diagonal');
    end
    
    a_name=sprintf('prec_Dinv');
    sd.add_var(a_name,diag(prec.Dinv), 'type', 'real');
    
    
    if (~isdiag(prec.E))
        
        error('prec.E is not diagonal');
    end
    
    a_name=sprintf('prec_E');
    sd.add_var(a_name,diag(prec.E), 'type', 'real');
    
    
end

sd.hl('extern REAL rho_inv;\n');
sd.cl(' REAL rho_inv = %f;\n', 1/rho);

sd.hl('extern REAL tau_inv;\n');
sd.cl(' REAL tau_inv = %f;\n', 1/tau);


%% adaptive step size

if(strcmp(settings.adaptive,'yes'))
    
    
    sd.hl('extern REAL sqrt_Mu;\n');
    sd.cl(' REAL sqrt_Mu = %f;\n', sqrt(Mu));
    
    sd.hl('extern REAL eps;\n');
    sd.cl(' REAL eps = %f;\n', eps);
    
    sd.hl('extern REAL svd_QT;\n');
    sd.cl(' REAL svd_QT = %f;\n', (1-max(eig(dat.Q))/2*max(svds(prec.T))));
    
    sd.hl('extern REAL prev_tau;\n');
    sd.cl(' REAL prev_tau = %f;\n', tau);
    
    sd.hl('extern REAL prev_rho;\n');
    sd.cl(' REAL prev_rho = %f;\n', rho);
    
    sd.hl('extern REAL prev_alpha;\n');
    sd.cl(' REAL prev_alpha = %f;\n', alpha);
    
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Previous Coder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  KKT solve step. Step 1 in Algorithm

%%%% First extract important matrix and also save variables

% Assumes that the KKT matrix is constant
% Q   = dat.Q;
% rho = rho;
% L   = dat.L;
% A   = dat.A;
%
% %%% in case no preconditioning, make scaling equals to 1
% if(strcmp(settings.precond,'no'))
%     settings.E = eye(size(L,1));
%     prec.E = settings.E;
% end
%
% if(strcmp(settings.precond,'yes'))
%     prec.E = settings.E;
% end
%
% %save mymat K;
%
% %%%% calculate and transform the problem in suite sparse format.
%
% if(strcmp(settings.adaptive,'yes'))
%
%
%     %%%% We wont use following variables so set them to zero
%     sd.define('nn_lp',0,'int');
%     sd.define('LNZ_ss', 0, 'int');
%     sd.define('ANZ_ss', 0, 'int');
%     sd.define('N_ss', 0, 'int');
%     A = 0;
%     sd.add_function(coderFunc_prefactor('custom_compute_prefactor',A,'none'));
%
%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%% copy paste starts from here
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     n = size(Q,1);
%     prox = prob.prox;
%     nProx = length(prox);
%     LL = zeros(n,n);
%     sum_l = zeros(n,1);
%     sum_l_base = zeros(n,1);
%     sum_L = [];
%     L = []; l = [];
%     nProxVars = zeros(nProx,1);
%     kk = 1;
%     for i = 1:nProx
%         proxInd{i}   = kk + [0:size(prox(i).L,1)-1];
%         ind          = proxInd{i};
%         sum_L        = [sum_L prox(i).L'*prec.E(ind,ind)'];
%         sum_l        = sum_l + rho*prox(i).L'*(prec.E(ind,ind)'*prec.E(ind,ind))*prox(i).l;
%         sum_l_base   = sum_l_base + prox(i).L'*(prec.E(ind,ind)'*prec.E(ind,ind))*prox(i).l;
%         l            = [l;prec.E(ind,ind)*prox(i).l]; % ATTENTION! L and l are left scaled
%         L            = [L;prec.E(ind,ind)*prox(i).L];
%         Q            = Q + rho*prox(i).L'*(prec.E(ind,ind)'*prec.E(ind,ind))*prox(i).L;
%         LL           = LL + (prec.E(ind,ind)*prox(i).L)'*(prec.E(ind,ind)*prox(i).L);
%         kk           = kk + size(prox(i).L,1);
%         nProxVars(i) = size(prox(i).L,1);
%     end
%     LL = full(LL);
%     dat.A = full(dat.A);
%     proxWeight = [prob.prox.weight]; % Weights for the prox functions
%
%
%     % Simultaneous diagonalization of KKT system
%
%     if (min(eig(LL))>0) && (isempty(dat.A))
%         H = dat.Q;   M = LL;
%         Lchol = chol(M,'lower');
%         C = Lchol \ (H*inv(Lchol'));
%         [V,D1] = qdwheig(C);
%         cas = 1;
%         sd.hl('#define adap_case_1\n');
%         sd.hl('extern const REAL rho1;\n');
%         sd.cl('const REAL rho1 = 1.0;\n');
%         sd.hl('extern const REAL rho2;\n');
%
%     elseif (min(eig(dat.Q))>0) && (isempty(dat.A))
%         H = LL;   M = dat.Q;
%         Lchol = chol(M,'lower');
%         C = Lchol \ (H*inv(Lchol'));
%         [V,D1] = qdwheig(C);
%         cas = 2;
%         sd.hl('#define adap_case_2\n');
%         sd.hl('extern const REAL rho1;\n');
%         sd.hl('extern const REAL rho2;\n');
%         sd.cl('const REAL rho2 = %1.0;\n');
%
%     else
%         error('Conditions for simultaneous diagonalization not satisfied. Adaptive ADMM not adviced - increased cost per iteration or try non-adaptive ADMM')
%         settings.adapt = 'no';
%
%     end
%     X = Lchol' \ V;
%
%
%     % Permute D - was not working otherwise
%     D1 = blkdiag(diag(nonzeros(diag(X'*H*X))),zeros(size(D1,1)-size(nonzeros(diag(X'*H*X)),1),size(D1,1)-size(nonzeros(diag(X'*H*X)),1)));
%     D1 = sparse(D1);
%
%
%     if(~isdiag(D1))
%         error('(D1 must be diagonal');
%     end
%
%
%
%
%     %%%%%% Add variable D1
%     D1_vec = full(diag(D1));
%     a_name=sprintf('D1_vec');
%     sd.add_var(a_name,D1_vec, 'type', 'real');
%
%     %%%%% Add variable X
% %     X_adapt_vec = vec(X');
% %     a_name=sprintf('X_adapt_vec');
% %     sd.add_var(a_name,X_adapt_vec, 'type', 'real');
%
%     %%%%% Add variable X'
%     %XTrans_adapt_vec = vec(X);
%     %a_name=sprintf('XTrans_adapt_vec');
%     %sd.add_var(a_name,XTrans_adapt_vec, 'type', 'real');
%
%
%     %%% add matrix vector multiplication
%     sd.add_function(...
%             coderFunc_times('custom_mult_XTrans', sparse(X'), 'Astr', 'XTrans_adapt_vec', 'method', Mat_Vec));
%     sd.add_function(...
%             coderFunc_times('custom_mult_X', sparse(X), 'Astr', 'X_adapt_vec', 'method', Mat_Vec));
%
%
%     if(strcmp(settings.precond,'yes'))
%
%         Ld = settings.E * dat.L;
%         sd.add_function(...
%             coderFunc_times('custom_mult_Ldtrans', sparse(Ld'), 'Astr', 'Ldtrans', 'method', Mat_Vec));
%
%
%     end
%
%
%     if(strcmp(settings.precond,'no'))
%      sd.add_function(...
%             coderFunc_times('custom_mult_Ltrans', sparse(dat.L'), 'Astr', 'Ltrans', 'method', Mat_Vec));
%     end
%
%
%
% end
%
% if(strcmp(settings.adaptive,'no'))
%
%     if(strcmp(settings.precond,'yes'))
%
%         Ld = settings.E * dat.L;
%         sd.add_function(...
%             coderFunc_times('custom_mult_Ldtrans', sparse(Ld'), 'Astr', 'Ldtrans', 'method', Mat_Vec));
%
%         K = [dat.Q+rho*Ld'*Ld A'; A zeros(size(A,1))];
%         mldivide = coderFunc_mldivide('custom_solve_kkt', K, 'Astr', 'KKT','method',Lin_Solve);
%         if strcmp(mldivide.desc,'ldl_ss')
%           sd.cl('double *Lx_ss;');
%           sd.cl('int *Li_ss;');
%           sd.hl('#define suitesparse_linsolve\n');
%
%         end
%
%  sd.add_function(mldivide);
%
%     end
%
%     if(strcmp(settings.precond,'no'))
%         sd.add_function(...
%             coderFunc_times('custom_mult_Ltrans', sparse(dat.L'), 'Astr', 'Ltrans', 'method', Mat_Vec));
%         K = [dat.Q+rho*L'*L A'; A zeros(size(A,1))];
%         mldivide = coderFunc_mldivide('custom_solve_kkt', K, 'Astr', 'KKT','method',Lin_Solve);
%         if strcmp(mldivide.desc,'ldl_ss')
%             sd.cl('double *Lx_ss;');
%             sd.cl('int *Li_ss;');
%             sd.hl('#define suitesparse_linsolve\n');
%
%         end
%
%         sd.add_function(mldivide);
%     end
%
%      LNZ=nnz(K)-nnz(triu(K));
% ANZ=nnz(triu(K));
% N=size(K,2);
%
% %%%% add varibles
%
% sd.define('nn_lp',size(K,2),'int');
% sd.define('LNZ_ss', LNZ, 'int');
% sd.define('ANZ_ss', ANZ, 'int');
% sd.define('N_ss', N, 'int');
%
%
% end
%
%
% %% Calculate prox function
%
% if(strcmp(settings.precond,'yes'))
%         %%% add ld variable here
%         ld = settings.E*dat.l;
%         a_name = sprintf('ld');
%         sd.add_var(a_name,ld,'type','real');
%     % workDual = L*x
%     sd.add_function(coderFunc_times('custom_mult_Ld', sparse(Ld), ...
%         'y_name', 'workDual', 'x_name', 'x', 'Astr', 'Ld', 'method', Mat_Vec));
%
%     % Evaluate prox functions y = prox(workDual)
%     sd.add_function(coderFunc_prox(dat));
%
%     %%%%%% Add variable E
%         E = settings.E;
%         E_vec = diag(E);
%         a_name=sprintf('E_vec');
%         sd.add_var(a_name,E_vec, 'type', 'real');
%
%         %%%%%% Add variable Einv
%         E_inv = inv(E);
%         Einv_vec = diag(E_inv);
%         a_name=sprintf('Einv_vec');
%         sd.add_var(a_name,Einv_vec, 'type', 'real');
%
% end
%
% if(strcmp(settings.precond,'no'))
%     % workDual = L*x
%     sd.add_function(coderFunc_times('custom_mult_L', sparse(dat.L), ...
%         'y_name', 'workDual', 'x_name', 'x', 'Astr', 'L', 'method', Mat_Vec));
%
%     % Evaluate prox functions y = prox(workDual)
%     sd.add_function(coderFunc_prox(dat));
%
% end
%
%
%
% %% for Residual calculation
%
% %%%%%% Add variable D1
%
%
% if(strcmp(settings.precond, 'yes'))
%
%     RES_mat = (Ld*settings.D_scale)';
%     sd.add_function(coderFunc_times('custom_mult_residual', RES_mat, ...
%         'y_name', 's', 'x_name', 'workDual', 'Astr', 'RES_mat', 'method', Mat_Vec));
%
%
% end

end

