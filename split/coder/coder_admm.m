% Generate custom c-code to solve ADMM
%% settings.E is such that yd = Ey where y is unscaled variables
%  settings.Dinv is such that and Dinv = inv(D) and Zp = Dz where z is unscaled variables
function sd = coder_admm(prob,settings)


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


if(strcmp(settings.precond, 'yes') && ~isfield(settings, 'E'))
    error('User must provide scaling matrix E for preconditioning')
end

if(strcmp(settings.precond, 'yes') && ~isfield(settings, 'D_scale'))
    error('User must provide scaling matrix D_scale for preconditioning')
end

if(isfield(settings, 'E') && ~isdiag(settings.E))
    error('Scaling matrix E must be diagonal')
end

if(isfield(settings, 'D_scale') && ~isdiag(settings.D_scale))
    error('Scaling matrix D_scale must be diagonal')
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

hard_mldivide = 0; % default is no PL or SOC. it is processors

FPGA_PL = 0;
if(isfield(settings, 'FPGA_PL'))
    
    if (settings.FPGA_PL == 1)
        FPGA_PL = 1;
        hard_mldivide = 2;
        
        
        if(~isfield(settings, 'Mat_Vec'))
            Mat_Vec = 'auto_FPGA';
        end
        
        if(~isfield(settings, 'Lin_Solve'))
            Lin_Solve = 'auto_FPGA';
            
        end
    end
    
end

FPGA_SoC = 0;
if(isfield(settings, 'FPGA_SoC'))
    
    if (settings.FPGA_SoC == 1)
        FPGA_SoC = 1;
        hard_mldivide = 1;
        
        
        if(~isfield(settings, 'Mat_Vec'))
            Mat_Vec = 'auto';
        end
        
        if(~isfield(settings, 'Lin_Solve'))
            Lin_Solve = 'auto_FPGA';
        end
        
    end
end

dat  = prob.coderData;
prox  = prob.prox;
sd   = splitData;


if FPGA_PL == 0
    
    sd.hl('#include "user_foo_data.h" \n');
    
else
    
    sd.hl('#include "foo_data.h" \n');
    
end


if(isfield(settings, 'primal_vars_x'))
    %% ifdef warm start x
    sd.hl('#define warm_start_x \n');
    a_name=sprintf('warm_x');
    sd.add_var(a_name,settings.primal_vars_x, 'type', 'real');
    
    
end

if(isfield(settings, 'primal_vars_y'))
    %% ifdef warm start y
    sd.hl('#define warm_start_y \n');
    a_name=sprintf('warm_y');
    if(strcmp(settings.precond,'yes'))
        sd.add_var(a_name,settings.E*settings.primal_vars_y, 'type', 'real');
    else
        sd.add_var(a_name,settings.primal_vars_y, 'type', 'real');
    end
end

if(isfield(settings, 'primal_vars_lam'))
    %% ifdef warm start lam
    sd.hl('#define warm_start_lam \n');
    a_name=sprintf('warm_lambda');
    sd.add_var(a_name,settings.primal_vars_lam, 'type', 'real');
end

if(isfield(settings, 'dual_vars_lam'))
    %% ifdef warm start lam
    sd.hl('#define warm_start_lambda \n');
    a_name=sprintf('warm_dual_lambda');
    sd.add_var(a_name,settings.dual_vars_lam, 'type', 'real');
end


if(strcmp(settings.adaptive,'yes'))
    
    gamma = (1 + sqrt(5))/2 - eps;
    sd.hl('#define alpha1\n');
    sd.define('gamma',gamma,'real');
else
    alpha = settings.relaxation;
    sd.hl('#define gamma1\n');
    sd.define('alpha',alpha,'real');
end

%%%%%% if preconditioned then add it in probdata.h file and also find
%%%%%% variables to save i.e. Ld, ld, D, D-1, E, E-1
if(strcmp(settings.precond,'yes'))
    sd.hl('#define precond\n');
end

if(strcmp(settings.adaptive,'yes'))
    sd.hl('#define adaptive\n');
end


if(isfield(settings, 'adaptive_everystep'))
    if(strcmp(settings.adaptive_everystep,'yes'))
        sd.hl('#define adaptive_everystep\n');
    end
end
sd.define('nParam',  size(dat.pL,2), 'int');
sd.define('nPrimal', size(dat.A,2),  'int');
sd.define('nDual',   size(dat.L,1),  'int');
sd.define('nEqCon',  size(dat.A,1),  'int');

if(~isfield(settings, 'rho'))
    rho = 1.0 ;
end

if(isfield(settings, 'rho'))
    rho = settings.rho;
    
end

if(strcmp(settings.adaptive,'yes') && isfield(settings, 'c_adapt_step'))
    
    if(isnumeric(settings.c_adapt_step))
        error('adaptive_step constant  must be numerical');
        
    else
        sd.hl('#define c_adapt_step %f \n',settings.c_adapt_step);
    end
    
else
    sd.hl('#define c_adapt_step %f \n',10);
end

rhoinv = 1/rho;

sd.hl('#define rho_init %f', rho);
sd.hl('#define rho_inv_init %f', rhoinv);


if(strcmp(settings.adaptive,'no'))
    sd.hl('extern  real rho;\n');
    %     sd.cl(' REAL rho = %f;\n', rho);
    sd.hl('extern  real rhoinv; \n');
    %     sd.cl(' REAL rhoinv = %f;', rhoinv);
end

if(strcmp(settings.adaptive,'yes'))
    sd.hl('extern real rho_tmp;\n');
    %     sd.cl(' REAL rho_tmp = %f;\n', rho);
    sd.hl('extern real rhoinv_tmp; \n');
    %     sd.cl(' REAL rhoinv_tmp = %f;', rhoinv);
end



%% custom_compute_parametric

% Compute: l = pL*par + l_const, etc

% Compute: l = pL*par + l_const, etc
if FPGA_PL == 1
    method_para = 'auto_FPGA';
    sd.add_function(coderFunc_parametric('custom_compute_parametric', dat, 'method', method_para ));
    
elseif FPGA_SoC == 1
    
    method_para = 'auto';
    sd.add_function(coderFunc_parametric('custom_compute_parametric', dat, 'method', method_para ));
    
else
    
    method_para = 'auto';
    sd.add_function(coderFunc_parametric('custom_compute_parametric', dat, 'method', method_para ));
    
end

%%  KKT solve step. Step 1 in Algorithm

if(~isfield(settings, 'latency'))
    settings.latency = 8;
    
end

if(~isfield(settings, 'paral'))
    settings.paral = 1;
end

%%%% First extract important matrix and also save variables

% Assumes that the KKT matrix is constant
Q   = dat.Q;
rho = rho;
L   = dat.L;
A   = dat.A;

%%% in case no preconditioning, make scaling equals to 1
if(strcmp(settings.precond,'no'))
    settings.E = eye(size(L,1));
    prec.E = settings.E;
end

if(strcmp(settings.precond,'yes'))
    prec.E = settings.E;
end


%%%% calculate and transform the problem in suite sparse format.

if(strcmp(settings.adaptive,'yes'))
    
    
    %%%% We wont use following variables so set them to zero
    sd.define('nn_lp',0,'int');
    sd.define('LNZ_ss', 0, 'int');
    sd.define('ANZ_ss', 0, 'int');
    sd.define('N_ss', 0, 'int');
    A = 0;
    sd.add_function(coderFunc_prefactor('custom_compute_prefactor',A,'none'));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% copy paste starts from here
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n = size(Q,1);
    prox = prob.prox;
    nProx = length(prox);
    LL = zeros(n,n);
    sum_l = zeros(n,1);
    sum_l_base = zeros(n,1);
    sum_L = [];
    L = []; l = [];
    nProxVars = zeros(nProx,1);
    kk = 1;
    for i = 1:nProx
        proxInd{i}   = kk + [0:size(prox(i).L,1)-1];
        ind          = proxInd{i};
        sum_L        = [sum_L prox(i).L'*prec.E(ind,ind)'];
        sum_l        = sum_l + rho*prox(i).L'*(prec.E(ind,ind)'*prec.E(ind,ind))*prox(i).l;
        sum_l_base   = sum_l_base + prox(i).L'*(prec.E(ind,ind)'*prec.E(ind,ind))*prox(i).l;
        l            = [l;prec.E(ind,ind)*prox(i).l]; % ATTENTION! L and l are left scaled
        L            = [L;prec.E(ind,ind)*prox(i).L];
        Q            = Q + rho*prox(i).L'*(prec.E(ind,ind)'*prec.E(ind,ind))*prox(i).L;
        LL           = LL + (prec.E(ind,ind)*prox(i).L)'*(prec.E(ind,ind)*prox(i).L);
        kk           = kk + size(prox(i).L,1);
        nProxVars(i) = size(prox(i).L,1);
    end
    
    
    LL = full(LL);
    dat.A = full(dat.A);
    proxWeight = [prob.prox.weight]; % Weights for the prox functions
    
    
    % Simultaneous diagonalization of KKT system
    
    if (min(eig(LL))>0) && (isempty(dat.A))
        H = dat.Q;   M = LL;
        Lchol = chol(M,'lower');
        C = Lchol \ (H*inv(Lchol'));
        [V,D1] = qdwheig(C);
        cas = 1;
        sd.hl('#define adap_case_1\n');
        sd.hl('extern const real rho1;\n');
        sd.cl('const REAL rho1 = 1.0;\n');
        sd.hl('extern const real rho2;\n');
        
    elseif (min(eig(dat.Q))>0) && (isempty(dat.A))
        H = LL;   M = dat.Q;
        Lchol = chol(M,'lower');
        C = Lchol \ (H*inv(Lchol'));
        [V,D1] = qdwheig(C);
        cas = 2;
        sd.hl('#define adap_case_2\n');
        sd.hl('extern const real rho1;\n');
        sd.hl('extern const real rho2;\n');
        sd.cl('const real rho2 = %1.0;\n');
        
    else
        error('Conditions for simultaneous diagonalization not satisfied. Adaptive ADMM not adviced - increased cost per iteration or try non-adaptive ADMM')
        settings.adapt = 'no';
        
    end
    X = Lchol' \ V;
    
    
    % Permute D - was not working otherwise
    D1 = blkdiag(diag(nonzeros(diag(X'*H*X))),zeros(size(D1,1)-size(nonzeros(diag(X'*H*X)),1),size(D1,1)-size(nonzeros(diag(X'*H*X)),1)));
    D1 = sparse(D1);
    
    
    if(~isdiag(D1))
        error('(D1 must be diagonal');
    end
    
    
    
    
    %%%%%% Add variable D1
    D1_vec = full(diag(D1));
    a_name=sprintf('D1_vec');
    sd.add_var(a_name,D1_vec, 'type', 'real');
    
    %%%%% Add variable X
    %     X_adapt_vec = vec(X');
    %     a_name=sprintf('X_adapt_vec');
    %     sd.add_var(a_name,X_adapt_vec, 'type', 'real');
    
    %%%%% Add variable X'
    %XTrans_adapt_vec = vec(X);
    %a_name=sprintf('XTrans_adapt_vec');
    %sd.add_var(a_name,XTrans_adapt_vec, 'type', 'real');
    
    
    %%% add matrix vector multiplication
    
    X_trans_name = coderFunc_times('custom_mult_XTrans', sparse(X'), 'Astr', 'XTrans_adapt_vec', 'method', Mat_Vec);
    sd.add_function(Xtrans_name);
    
    if strcmp(X_trans_name.desc, 'FPGA_matvec_sparse')
        
        sd.hl('#include "user_custom_mult_Xtrans_sparse_mv_mult.h" \n');
        
        
    end
    
    X_name = coderFunc_times('custom_mult_X', sparse(X), 'Astr', 'X_adapt_vec', 'method', Mat_Vec);
    sd.add_function(X_name);
    
    if strcmp(X_name.desc, 'FPGA_matvec_sparse')
        
        sd.hl('#include "user_custom_mult_X_sparse_mv_mult.h" \n');
        
        
    end
    
    if(strcmp(settings.precond,'yes'))
        
        Ld = settings.E * dat.L;
        
        L_trans_name = coderFunc_times('custom_mult_Ldtrans', sparse(Ld'), 'Astr', 'Ldtrans', 'method', Mat_Vec);
        sd.add_function(L_trans_name);
        
        if strcmp(L_trans_name.desc, 'FPGA_matvec_sparse')
            
            sd.hl('#include "user_custom_mult_Ldtrans_sparse_mv_mult.h" \n');
            
            
        end
        
    end
    
    
    if(strcmp(settings.precond,'no'))
        L_trans_name = coderFunc_times('custom_mult_Ltrans', sparse(dat.L'), 'Astr', 'Ltrans', 'method', Mat_Vec);
        sd.add_function(L_trans_name);
        
        if strcmp(L_trans_name.desc, 'FPGA_matvec_sparse')
            
            sd.hl('#include "user_custom_mult_Ltrans_sparse_mv_mult.h" \n');
            
            
        end
    end
    
    
    
end

if(strcmp(settings.adaptive,'no'))
    
    if(strcmp(settings.precond,'yes'))
        
        %% printing
        disp('proxil is ')
        prox(1).l
        %%
        Ld = settings.E * dat.L;
        
        L_trans_name = coderFunc_times('custom_mult_Ldtrans', sparse(Ld'), 'Astr', 'Ldtrans', 'method', Mat_Vec);
        sd.add_function(L_trans_name);
        
        if strcmp(L_trans_name.desc, 'FPGA_matvec_sparse')
            
            sd.hl('#include "user_custom_mult_Ldtrans_sparse_mv_mult.h" \n');
            
            
        end
        
        K = [dat.Q+rho*Ld'*Ld A'; A zeros(size(A,1))];
        spy(K)
        
        if(~isfield(settings, 'latency'))
            latency = 8;
        end
        
        if(~isfield(settings, 'paral'))
            paral = floor(size(K,1)/3);
        end
        
        nPrimal_solve = size(dat.Q,1);
        mldivide = coderFunc_mldivide('custom_solve_kkt', K, 'Astr', 'KKT','method',Lin_Solve, 'nPrimal', nPrimal_solve, 'lat' ,settings.latency, 'paral', settings.paral, 'hard', hard_mldivide);
        
        
        if strcmp(mldivide.desc,'ldl_ss')
            sd.cl('double *Lx_ss;');
            sd.cl('int *Li_ss;');
            sd.hl('#define suitesparse_linsolve\n');
            
            
        end
        
        if strcmp(mldivide.desc,'ldl_lp')
            sd.cl('__CLPK_integer ipiv[nn_lp];');
            
            sd.hl('#define lapack_linsolve\n');
            
        end
        
        
        
        sd.add_function(mldivide);
        
    end
    
    if(strcmp(settings.precond,'no'))
        
        L_trans_name = coderFunc_times('custom_mult_Ltrans', sparse(dat.L'), 'Astr', 'Ltrans', 'method', Mat_Vec);
        sd.add_function(L_trans_name);
        
        if strcmp(L_trans_name.desc, 'FPGA_matvec_sparse')
            
            sd.hl('#include "user_custom_mult_Ltrans_sparse_mv_mult.h" \n');
            
            
        end
        
        K = [dat.Q+rho*L'*L A'; A zeros(size(A,1))];
        
        nPrimal_solve = size(dat.Q,1);
        
        if(~isfield(settings, 'latency'))
            latency = 8;
        end
        
        if(~isfield(settings, 'paral'))
            paral = floor(size(K,1)/3);
        end
        
        mldivide = coderFunc_mldivide('custom_solve_kkt', K, 'Astr', 'KKT','method',Lin_Solve, 'nPrimal', nPrimal_solve, 'lat' , settings.latency, 'paral', settings.paral, 'hard', hard_mldivide);
        
        if strcmp(mldivide.desc,'ldl_ss')
            sd.cl('double *Lx_ss;');
            sd.cl('int *Li_ss;');
            sd.hl('#define suitesparse_linsolve\n');
            
        end
        
        if strcmp(mldivide.desc,'ldl_lp')
            sd.cl('__CLPK_integer ipiv[nn_lp];');
            
            sd.hl('#define lapack_linsolve\n');
            
        end
        
        
        
        sd.add_function(mldivide);
        
    end
    
    if (strcmp(mldivide.desc,'invert_FPGA_tree') || strcmp(mldivide.desc,'invert_FPGA_MAC'))
        
        sd.hl('#include "user_mv_mult.h" \n');
    end
    
    
    
    LNZ=nnz(K)-nnz(triu(K));
    ANZ=nnz(triu(K));
    N=size(K,2);
    
    %%%% add varibles
    
    sd.define('nn_lp',size(K,2),'int');
    sd.define('LNZ_ss', LNZ, 'int');
    sd.define('ANZ_ss', ANZ, 'int');
    sd.define('N_ss', N, 'int');
    
    
end



%% Calculate prox function

sd.hl('extern real prox_var; \n');
sd.cl(' real prox_var = %f;', 0.0);

if(strcmp(settings.precond,'yes'))
    
    
    ld = settings.E*prox.l;
    prox.l
    %%%%%%%%% copy paste ends
    %%% add ld variable here
    
    a_name = sprintf('ld');
    sd.add_var(a_name,ld,'type','real');
    % workDual = L*x
    
    
    L_times_name = (coderFunc_times('custom_mult_Ld', sparse(Ld), ...
        'y_name', 'workDual', 'x_name', 'x', 'Astr', 'Ld', 'method', Mat_Vec));
    sd.add_function = (L_times_name);
    % Evaluate prox functions y = prox(workDual)
    sd.add_function(coderFunc_prox(dat));
    
    %%%%%% Add variable E
    E = settings.E;
    E_vec = diag(E);
    a_name=sprintf('E_vec');
    sd.add_var(a_name,E_vec, 'type', 'real');
    
    %%%%%% Add variable Einv
    E_inv = inv(E);
    Einv_vec = diag(E_inv);
    a_name=sprintf('Einv_vec');
    sd.add_var(a_name,Einv_vec, 'type', 'real');
    
end

if(strcmp(settings.precond,'no'))
    % workDual = L*x
    
    L_times_name =(coderFunc_times('custom_mult_L', sparse(dat.L), ...
        'y_name', 'workDual', 'x_name', 'x', 'Astr', 'L', 'method', Mat_Vec));
    sd.add_function(L_times_name);
    % Evaluate prox functions y = prox(workDual)
    sd.add_function(coderFunc_prox(dat));
    
end

if strcmp(L_times_name.desc, 'FPGA_matvec_sparse')
    
    sd.hl('#include "user_custom_mult_L_sparse_mv_mult.h" \n');
    
    
end



%% for Residual calculation

%%%%%% Add variable D1


if(strcmp(settings.precond, 'yes'))
    
    RES_mat = (Ld*settings.D_scale)';
    sd.add_function(coderFunc_times('custom_mult_residual', RES_mat, ...
        'y_name', 's', 'x_name', 'workDual', 'Astr', 'RES_mat', 'method', Mat_Vec));
    
    
end

end

