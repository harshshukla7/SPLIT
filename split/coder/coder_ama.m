% Generate custom c-code to solve AMA
function sd = coder_ama(prob,settings)
%% settings.E is such that yd = Ey where y is unscaled variables
%  settings.Dinv is such that and Dinv = inv(D_scale) and Zp = D_scale*z where z is unscaled variables
%  Matrix-vector product options :  'forloops', 'sparse', 'blas', 'ss', exhaustive_gen'
%  Linear solve options          :  'invert','ldl','ldl_lp','ldl_ss'


if(~isfield(settings, 'precond'))
    settings.precond = 'no';
end



if(strcmp(settings.precond,'yes') && strcmp(settings.precond,'no'))
    error('precond value must be either "yes" or "no"');
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

if(strcmp(settings.precond, 'yes'))
    
    prec.E = settings.E;
    prec.Dinv = settings.D_scale;
    warning('Appropariate Dimension of the scaling matrix is not checked by SPLIT. Make sure its of correct dimension.')
end

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
    Lin_Solve = 'ldl_ss';
end

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
            Lin_Solve = 'auto_SoC';
        end
        
    end
end


dat  = prob.coderData;
sd   = splitData;
prox  = prob.prox;

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
    sd.add_var(a_name,settings.primal_vars_y, 'type', 'real');
end

if(isfield(settings, 'primal_vars_lam'))
    %% ifdef warm start lam
    sd.hl('#define warm_start_lam \n');
    a_name=sprintf('warm_lambda');
    sd.add_var(a_name,settings.primal_vars_lam, 'type', 'real');
end


%%%%%% if preconditioned then add it in probdata.h file and also find
%%%%%% variables to save i.e. Ld, ld, D, D-1, E, E-1
if(strcmp(settings.precond,'yes'))
    sd.hl('#define precond\n');
end

if(isfield(settings, 'restart'))
    if(strcmp(settings.restart, 'yes'))
        sd.hl('#define adaptive_restart\n');
    end
    
end

%%%% calculate and transform the problem in suite sparse format.

sd.define('nParam',  size(dat.pL,2), 'int');
sd.define('nPrimal', size(dat.A,2),  'int');
sd.define('nDual',   size(dat.L,1),  'int');
sd.define('nEqCon',  size(dat.A,1),  'int');

Q   = dat.Q;
L   = dat.L;
A   = dat.A;

% Test: tune the stepsize rho:

if strcmp(settings.precond, 'yes')
    sigma_f = min(eig(prec.Dinv'*dat.Q*prec.Dinv));
    
else
    sigma_f = min(eig(dat.Q));
end


if sigma_f == 0
    error('Cannot use FAMA - objective not strongly convex') % XXX to be checked again; condition on nullspace sufficient
end

eigenvalue_LL = max(eig(L'*L));


if(~isfield(settings, 'rho'))
    rho = sigma_f/eigenvalue_LL(end);
end


if(isfield(settings, 'rho'))
    rho = settings.rho;
end


rhoinv = 1/rho;


% sd.hl('\n const real rho = %f;', rho);
% sd.hl('const real rhoinv = %f;', rhoinv);
% sd.hl(' real prox_var = %f;\n', 0.0);

sd.hl('#define rho_init %f', rho);
sd.hl('#define rho_inv_init %f', rhoinv);
sd.hl('#define prox_var_init %f', 0.0);



sd.hl('extern real rho;\n');
%  sd.cl('const REAL rho = %f;\n', rho);
sd.hl('extern real rhoinv; \n');
%  sd.cl('const REAL rhoinv = %f;', rhoinv);
sd.hl('extern real prox_var; \n');
%  sd.cl(' REAL prox_var = %f;', 0.0);


%% custom_compute_parametric

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
%% kktRHS[1:nPrimal] = L'*workDual
% y = L'*x

if(strcmp(settings.precond,'yes'))
    
    Ld = settings.E * dat.L;
    
    Lt_times_name = coderFunc_times('custom_mult_Ldtrans', (Ld'), 'Astr', 'Ldtrans', 'method', Mat_Vec);
    sd.add_function(Lt_times_name);
        
    
    
end

if(strcmp(settings.precond,'no'))
    
     Lt_times_name = coderFunc_times('custom_mult_Ltrans', (dat.L'), 'Astr', 'Ltrans', 'method', Mat_Vec);
    
    sd.add_function(Lt_times_name);
    
end

if strcmp(Lt_times_name.desc, 'FPGA_matvec_sparse')
   
        sd.hl('#include "user_custom_mult_Ltrans_sparse_mv_mult.h" \n');
        
    
end
%%  checking matrix vector multpilication
%
% TP1 = [1,2 ; 3, 4];
%
%
% sd.add_function(...
%             coderFunc_times('custom_mult_TP1', TP1, 'Astr', 'TP1'));

%% Solve the KKT system

K = [Q A'; A zeros(size(A,1))]; %% KKT system matrix


if(~isfield(settings, 'latency'))
    settings.latency = 8;
    
end

if(~isfield(settings, 'paral'))
    settings.paral = 1;
end


LNZ=nnz(K)-nnz(triu(K));
ANZ=nnz(triu(K));
N=size(K,2);
nPrimal_solve = size(Q,1);
%%%% add the varibles

sd.define('nn_lp',size(K,2),'int');
sd.define('LNZ_ss', LNZ, 'int');
sd.define('ANZ_ss', ANZ, 'int');
sd.define('N_ss', N, 'int');

%%%%% edited part ends here.
% if strcmp(Lin_Solve, 'ldl_ss')
% sd.cl('double *Lx_ss;');
% sd.cl('int *Li_ss;');
% end

%sd.add_function(coderFunc_mldivide('custom_solve_kkt', K, 'Astr', 'KKT','method',Lin_Solve));

mldivide = coderFunc_mldivide('custom_solve_kkt', K, 'Astr', 'KKT','method', Lin_Solve, 'nPrimal', nPrimal_solve, 'lat' , settings.latency, 'paral', settings.paral, 'hard', hard_mldivide );


if strcmp(mldivide.desc,'ldl_ss')
    sd.cl('double *Lx_ss;');
    sd.cl('int *Li_ss;');
    sd.hl('#define suitesparse_linsolve\n');
    
end

if strcmp(mldivide.desc,'ldl_lp')
    sd.cl('__CLPK_integer ipiv[nn_lp];');
    
    sd.hl('#define lapack_linsolve\n');
    
end

if ((strcmp(mldivide.desc,'invert_FPGA_tree') || strcmp(mldivide.desc,'invert_FPGA_MAC')) && (FPGA_PL == 1))
    
    sd.hl('#include "user_mv_mult.h" \n');
end

sd.add_function(mldivide);


%%
if(strcmp(settings.precond,'yes'))
    %%% add ld variable here
    ld = settings.E*dat.l;
    a_name = sprintf('ld');
    sd.add_var(a_name,ld,'type','real');
    % workDual = L*x
    L_times_name = coderFunc_times('custom_mult_Ld', (Ld), ...
        'y_name', 'workDual', 'x_name', 'x', 'Astr', 'Ld', 'method', Mat_Vec);
    sd.add_function(L_times_name);
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
    L_times_name = coderFunc_times('custom_mult_L', full(dat.L), ...
        'y_name', 'workDual', 'x_name', 'x', 'Astr', 'L','method',Mat_Vec, 'loop_reset', 0);
    sd.add_function(L_times_name);
    
end


if strcmp(L_times_name.desc, 'FPGA_matvec_sparse')
   
        sd.hl('#include "user_custom_mult_L_sparse_mv_mult.h" \n');
        
    
end

%% Evaluate prox functions y = prox(workDual)
sd.add_function(coderFunc_prox(dat));


%% For residual calculation
if(strcmp(settings.precond, 'yes'))
    
    RES_mat = (Ld*settings.D_scale)';
    sd.add_function(coderFunc_times('custom_mult_residual', RES_mat, ...
        'y_name', 's', 'x_name', 'workDual', 'Astr', 'RES_mat', 'method', Mat_Vec));
    
    
end




end

