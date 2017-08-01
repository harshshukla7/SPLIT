%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast gradient method and gradient method
% Example: ball on table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

%% Ball on table
% x0 = [10;10;10;10];
% N = 10;
% dat = ball_on_table(x0,N);

%% Quadrotor
N = 15;
x0 = 0.5*[-1,(10/180)*pi,-(10/180)*pi,(120/180)*pi,0,0,0]';
x0 = [-1,(10/180)*pi,-(10/180)*pi,(120/180)*pi,0,0,0]';
dat = quadrotor_model(x0,N);


%% Build matrices for fast gradient method
A_g = [];
B_g = zeros((dat.N)*size(dat.B,1),(dat.N-1)*size(dat.B,2));
Q_g = kron(dat.Q,eye(dat.N));
R_g = kron(dat.R,eye(dat.N-1));

B_help = [];

for j = 1:dat.N-1
        B_help = [dat.A^(j-1)*dat.B, B_help];
end

for i = 1: dat.N 
    A_g = [A_g' dat.A^(i-1)']';
    
    if i > 1
        B_g(1+(i-1)*size(dat.B,1):(i)*size(dat.B,1),1:(i-1)*size(dat.B,2)) = B_help(:,end-(i-1)*size(dat.B,2)+1:end);
    end
end

H_g = B_g'*Q_g*B_g+R_g;
M1_g = A_g'*Q_g*B_g;
M2_g = A_g'*Q_g*A_g;

dat.H = H_g;
dat.h = dat.x0'*M1_g;

%% parameter

dat.mu = min(eig(H_g));
dat.L = max(eig(H_g));
dat.alpha0 = (1-sqrt(dat.mu/dat.L))/5*4;
dat.gamma0 = dat.alpha0*(dat.alpha0*dat.L-dat.mu)/(1-dat.alpha0);

%% Yalmip solution
sol_u = yal_solver_for_FGM(dat);
u_star = sol_u.u(:);
obj_star = sol_u.obj;

%% FGM algorithm

step = 100;
alpha = dat.alpha0;
alpha_1 = 0;
beta = 0;
u_fg = zeros((dat.N-1)*dat.m,step);
u_fg_hat = zeros((dat.N-1)*dat.m,step);
u_min = dat.umin;
u_max = dat.umax;
diff_u_fg = zeros(step,1);
diff_obj_fg = zeros(step,1);
diff_u_fg(1) = (u_fg(:,1)-u_star)'*(u_fg(:,1)-u_star);
diff_obj_fg(1) = 1/2*u_fg(:,1)'*H_g*u_fg(:,1) + dat.x0'*M1_g*u_fg(:,1) + 1/2*dat.x0'*M2_g*dat.x0-obj_star;


%% Yalmip solution
u_test = sdpvar(dat.m*(dat.N-1),1);

con_test = [];
con_test = con_test + set(u_min <= u_test <= u_max);

obj_test = 1/2*u_test'*H_g*u_test + dat.x0'*M1_g*u_test;

sol2 = solvesdp(con_test,obj_test);
sol_Dense.u = double(u_test);

u_diff_stant_dense = (u_star-sol_Dense.u)'*(u_star-sol_Dense.u);

%% FGD
for i = 2:step
    
    % step 1
    u_fg(:,i) = u_fg_hat(:,i-1)-1/dat.L*(H_g*u_fg_hat(:,i-1)+M1_g'*dat.x0);
    
    % step 2
    u_fg(:,i) = max(u_min,u_fg(:,i));
    u_fg(:,i) = min(u_max,u_fg(:,i));
    
    % save for plot
    diff_u_fg(i) = (u_fg(:,i)-sol_Dense.u)'*(u_fg(:,i)-sol_Dense.u);
    diff_obj_fg(i) = 1/2*u_fg(:,i)'*H_g*u_fg(:,i) + dat.x0'*M1_g*u_fg(:,i) + 1/2*dat.x0'*M2_g*dat.x0-obj_star;
 
    % step 3 : acceleratiom step
    poly_alpha = [1;(alpha^2-dat.mu/dat.L);-alpha^2];
    roots_alpha = roots(poly_alpha);
    
    if roots_alpha(1) < 1 && roots_alpha(1) > 0
        alpha_1 = roots_alpha(1);
    else 
        alpha_1 = roots_alpha(2);
    end
    
    beta = alpha*(1-alpha)/(alpha^2+alpha_1);
    
    u_fg_hat(:,i) = u_fg(:,i) + beta*(u_fg(:,i) - u_fg(:,i-1));
    
    alpha = alpha_1;

end



%% Gradient methods

u_g = zeros((dat.N-1)*dat.m,step);
u_min = -ones((dat.N-1)*dat.m,1);
u_max = ones((dat.N-1)*dat.m,1);
diff_u_g = zeros(step,1);
diff_obj_g = zeros(step,1);
diff_u_g(1) = (u_g(:,1)-u_star)'*(u_g(:,1)-u_star);
diff_obj_g(1) = 1/2*u_g(:,1)'*H_g*u_g(:,1) + dat.x0'*M1_g*u_g(:,1) + 1/2*dat.x0'*M2_g*dat.x0-obj_star;

for i = 2:step
    
    % step 1
    u_g(:,i) = u_g(:,i-1)-1/dat.L*(H_g*u_g(:,i-1)+M1_g'*dat.x0);
    % step 2
    u_g(:,i) = max(u_min,u_g(:,i));
    u_g(:,i) = min(u_max,u_g(:,i));
    % save
    diff_u_g(i) = (u_g(:,i)-sol_Dense.u)'*(u_g(:,i)-sol_Dense.u);
    diff_obj_g(i) = 1/2*u_g(:,i)'*H_g*u_g(:,i) + dat.x0'*M1_g*u_g(:,i) + 1/2*dat.x0'*M2_g*dat.x0-obj_star;
 
end


% %% Pre-conditioning
% C_g = chol(H_g);
% 
% dat.kapa_old = max(eig(C_g'*C_g))/min(eig(H_g));
% 
% E = sdpvar(size(H_g,1));
% r = sdpvar(1,1,'full');
% 
% con = [];
% obj = r;
% 
% con = con + set(H_g >= dat.mu*E >= 0);
% con = con + set([E C_g; C_g' r*eye(size(C_g,1))] >= 0);
% 
% opt = sdpsettings('solver','sdpt3');
% %opt = sdpsettings('solver','sedumi'); 
% solvesdp(con,obj,opt);
% 
% Sol.E = double(E);
% Sol.r = double(r);
% 
% P = inv(chol(Sol.E));
% 
% dat.kapa_new = max(eig(P'*H_g*P))/min(eig(P'*H_g*P));
% 
% %% Fast gradient with pre-conditioning
% % u = Pu_p, H_g_p = P'*H_g*P, u_min_p = inv(P)*u_min, u_max_p = inv(P)*u_max, M1_g_p =
% % M1_g*P, u_star = P*u_star_p
% 
% H_g_p = P'*H_g*P;
% dat.H_g_p = H_g_p;
% u_min_p = inv(P)*u_min;
% u_max_p = inv(P)*u_max;
% dat.u_min_p = u_min_p;
% dat.u_max_p = u_max_p;
% u_star_p = inv(P)*u_star;
% dat.mu_p = min(eig(dat.H_g_p));
% dat.L_p = max(eig(dat.H_g_p));
% dat.alpha0_p = (1-sqrt(dat.mu_p/dat.L_p))/5*4;
% dat.gamma0_p = dat.alpha0_p*(dat.alpha0_p*dat.L_p-dat.mu_p)/(1-dat.alpha0_p);
% M1_g_p = M1_g*P;
% dat.x0_p = dat.x0;
% 
% %%
% u_yal_p = sdpvar(size(u_star,1),size(u_star,2));
% 
% con_p = [];
% con_p = con_p + set(u_min <= P*u_yal_p <= u_max);
% 
% obj_p = 1/2*u_yal_p'*H_g_p*u_yal_p + dat.x0_p'*M1_g_p*u_yal_p;
% 
% solvesdp(con_p,obj_p);
% 
% 
% %%
% alpha_p = dat.alpha0_p;
% alpha_p_1 = 0;
% beta_p = 0;
% u_fg_p = zeros((dat.N-1)*dat.m,step);
% u_fg_hat_p = zeros((dat.N-1)*dat.m,step);
% diff_u_fg_p = zeros(step,1);
% diff_obj_fg_p = zeros(step,1);
% diff_u_fg_p(1) = (u_fg_p(:,1)-u_star_p)'*(u_fg_p(:,1)-u_star_p);
% diff_obj_fg_p(1) = 1/2*u_fg_p(:,1)'*H_g_p*u_fg_p(:,1) + dat.x0_p'*M1_g_p*u_fg_p(:,1) + 1/2*dat.x0_p'*M2_g*dat.x0_p-obj_star;
% 
% for i = 2:step
%     
%     u_fg_p(:,i) = u_fg_hat_p(:,i-1)-1/dat.L_p*(H_g_p*u_fg_hat_p(:,i-1)+M1_g_p'*dat.x0_p);
%     
% %     temp_p = max(u_min,P*u_fg_p(:,i));
% %     temp_p = min(u_max,temp_p);
%     temp_p = min([u_max;-u_min],[P;-P]*u_fg_p(:,i));
%     u_fg_p(:,i) = inv(P)*temp_p;
%     
%     diff_u_fg_p(i) = (u_fg_p(:,i)-u_star_p)'*(u_fg_p(:,i)-u_star_p);
%     diff_obj_fg_p(i) = 1/2*u_fg_p(:,i)'*H_g_p*u_fg_p(:,i) + dat.x0_p'*M1_g_p*u_fg_p(:,i) + 1/2*dat.x0_p'*M2_g*dat.x0_p-obj_star;
%  
%     poly_alpha_p = [1;(alpha_p^2-dat.mu_p/dat.L_p);-alpha_p^2];
%     roots_alpha_p = roots(poly_alpha_p);
%     
%     if roots_alpha_p(1) < 1 && roots_alpha_p(1) > 0
%         alpha_p_1 = roots_alpha_p(1);
%     else 
%         alpha_p_1 = roots_alpha_p(2);
%     end
%     
%     beta_p = alpha_p*(1-alpha_p)/(alpha_p^2+alpha_p_1);
%     
%     u_fg_hat_p(:,i) = u_fg_p(:,i) + beta_p*(u_fg_p(:,i) - u_fg_p(:,i-1));
%     
%     alpha_p = alpha_p_1;
% 
% end

%% plot
figure(1);grid on;

leg1(1) = semilogy([1:step],diff_u_g,'Color','Blue','LineWidth',1.5);hold on;
leg1(2) = semilogy([1:step],diff_u_fg,'Color','green','LineWidth',1.5);grid on;
%leg1(3) = semilogy([1:step],diff_u_fg_p,'Color','red','LineWidth',1.5);grid on;
figleg1 = legend(leg1,'|uk-u*| of gradient method','|uk-u*| of fast gradient method','|vk-v*| of fast gradient method with preconditionning');















































%% Quadsrotor
% dat.n = 7;
% dat.m = 4;
% dat.N = 8;
% 
% dat.A =  [1.0000         0         0         0         0         0         0;
%                0    1.0000         0         0    0.1000         0         0;
%                0         0    1.0000         0         0    0.1000         0;
%                0         0         0    1.0000         0         0    0.1000;
%                0         0         0         0    1.0000         0         0;
%                0         0         0         0         0    1.0000         0;
%                0         0         0         0         0         0    1.0000];
%            
% dat.B = [0.3500    0.3500    0.3500    0.3500;
%               0    0.0028         0   -0.0028;
%         -0.0028         0    0.0028         0;
%          0.0037   -0.0037    0.0037   -0.0037;
%               0    0.0560         0   -0.0560;
%         -0.0560         0    0.0560         0;
%          0.0733   -0.0733    0.0733   -0.0733];
% 
% us = [0.7007;
%       0.7007;
%       0.7007;
%       0.7007];
% 
% % Constraints
% dZMax = 1;
% angleMax =  [1;1]*(10/180)*pi;
% dangleMax =  [1;1]*(15/180)*pi;
% dgammaMax =  1*(60/180)*pi;
% 
% x1Max1 = [dZMax;angleMax];
% x1Max2 = [dangleMax;dgammaMax];
% 
% u1Min = zeros(dat.m,1)-us;
% u1Max = ones(dat.m,1)-us;
% 
% dat.C = zeros(dat.n*2+dat.m*2,dat.n);
% dat.C(1:dat.n,1:dat.n) = eye(dat.n);
% dat.C(dat.n+1:dat.n*2,1:dat.n) = -eye(dat.n);
% 
% dat.D = zeros(dat.n*2+dat.m*2,dat.m);
% dat.D(dat.n*2+1:dat.n*2+dat.m,1:dat.m) = eye(dat.m);
% dat.D(dat.n*2+dat.m+1:dat.n*2+2*dat.m,1:dat.m) = -eye(dat.m);
% 
% dat.c = [x1Max1;100;x1Max2;x1Max1;100;x1Max2;u1Max;u1Max];
%  
% %Stage cost
% dat.Q =diag([10, 60, 20, 2,1,1,1]);%diag([100, 6000, 2000, 2,0.1,0.1,0.1]);% diag([10,1000,1000,100,0.1,0.1,0.1]);
% dat.R = diag([1,1,1,1]);

%dat.x0 = [-1,(10/180)*pi,-0,(120/180)*pi,0,0,0]';
