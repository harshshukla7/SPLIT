function u_fg = solve_FGM(FGM, dat)
%% FGM algorithm
% ----------------------------------------------------------------------- %
%  TODO: unify the nomenclature, i.e. symbols in the pseudocode are not
%  identical with symbols in Ye's formulation
%  algorithm:
%       step1: projection
%              z_K = proj(y - 1/L*gradf(y))
%
%       step2: compute alpha such that \alpha_i \ in (0, 1) according to the
%       following equation: 
%              alpha_i+1^2 = (1 - alpha_i+1) alpha_i^2 + mu alpha_i+1/L
%
%       step3: compute beta_i, such that:
%              alpha_i(1-alpha_i)/alpha_i^2+alpha_i+1
%
%       step4: update y_i+1 = z_i+1 + beta_i(z_i+1-z_i)          
% ----------------------------------------------------------------------- %
% 

% initialize parameters used in the computation 
alpha_i = FGM.alpha0;
alpha_i_1 = 0;
beta_i = 0;
  
% memory allocations
u_fg = zeros((dat.N-1)*dat.m, 1);
u_fg_prev = zeros((dat.N-1)*dat.m, 1);
u_fg_hat = zeros((dat.N-1)*dat.m, 1);
u_fg_hat_prev = zeros((dat.N-1)*dat.m, 1);
    
% box constraints
u_min = -ones((dat.N-1)*dat.m, 1);
u_max = ones((dat.N-1)*dat.m, 1);

for i = 1:dat.step
    
    % projection
    u_fg = u_fg_hat_prev - 1/FGM.L * (FGM.H_g*u_fg_hat_prev + FGM.M1_g' * dat.x0);
    u_fg = max(u_min, u_fg);
    u_fg = min(u_max, u_fg);
    
    % compute the roots - in C it should be good, since it is just a
    % a quadratic equation
    poly_alpha = [1; (alpha_i^2-FGM.mu/FGM.L);-alpha_i^2];
    roots_alpha = roots(poly_alpha);
    
    % I need
    if roots_alpha(1) < 1 && roots_alpha(1) > 0
       alpha_i_1 = roots_alpha(1);
    else 
       alpha_i_1 = roots_alpha(2);
    end
    
    % compute
    beta_i = alpha_i*(1-alpha_i)/(alpha_i^2+alpha_i_1);
    
    % update of something - I have no idea, though what is it
    u_fg_hat = u_fg + beta_i*(u_fg - u_fg_prev);
   
    % save previous values
    alpha_i = alpha_i_1;
    u_fg_prev = u_fg;
    u_fg_hat_prev = u_fg_hat;
        
end
    