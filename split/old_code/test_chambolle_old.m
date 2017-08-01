function [Z,z,jjj] = test_chambolle(NO)

%% Solves the problem  minimize .5 * z'*H*z + l(T*z+t)
%%                     subject to A*z = b,
%% where l(.) is one of the functions presented below, according to the case chosen (NO)
randn('state',0);
rand('state',0);

MAXITER = 1e5;

%% Construct main part of the problem
n = 3;
A = randn(n,2*n); b = randn(n,1);
H1 = randn(2*n,2*n);
H = H1'*H1;
sigma_G = min(eig(H));
     
%% Run CPI 
switch NO
    
    case 1 %% l(T*z+t) = \delta+(T*z+t)
        disp('Nonegative orthant')
        %% Initialiations
        T = randn(3*n,2*n);   t = .18 - randn(3*n,1);
        K = [T; A];
        k = [t; b];
        no.con.all = size(K,1);
        Mu = max(svds(K)); % Bound on K operator
        theta = 1; sigma = 1*.9/Mu; tau = (1/1)*sigma;
        p = zeros(3*n,1); old.p = p; tilde.p = p; % has the size of T
        nu = zeros(n,1); old.nu = nu; tilde.nu = nu;
        z = zeros(2*n,1);   old.z = z;  bar.z = z;
        
        % precompute some stuff man
        matrix_inverse = inv(H+eye(size(H,1))*(1/tau));
        
        %% Main
        for jjj = 1:MAXITER 
            %% step 1
            old.p = p;
            old.nu = nu;
            tilde.p = p + sigma*(T*bar.z+t);
            nu = nu + sigma*(A*bar.z-b);

            p = proj_negative(tilde.p);

            %% step 2
            old.z = z;
            z = matrix_inverse * ((1/tau)*z-(K'*[p;nu]));

            %% step 3
            bar.z = z + theta*(z-old.z);

%             res.s  = tau*theta \ ( z-old.z);
%             res.r  = theta*K*(bar.z-old.z) + sigma \ ([p;nu]-[old.p;old.nu]);
            res.s  = 1/tau*theta*( z-old.z);
            res.r  = theta*K*(bar.z-old.z) + 1/sigma*([p;nu]-[old.p;old.nu]);
            
            rDual(jjj) = norm(res.s);    rPrimal(jjj) = norm(res.r);
            if (rDual(jjj) <= 1e-5 && rPrimal(jjj) <= 1e-5)
                    break;
            end
        end
        %%--------CVX---------%%
        cvx_begin
        cvx_solver 'sedumi'
            variables Z(2*n,1);
            obj = .5*Z'*H*Z;
            minimize ( obj )
            subject to 
            A*Z == b;
            T*Z+t >= 0;
        cvx_end
        fprintf('Chambolle solution\n'); disp(z)
        fprintf('CVX solution\n'); disp(Z)
        fprintf('Inequality constraints active:\t');
        if (any(T*Z+t)<=1e-5) 
            disp('YES') 
        else disp('NO') 
        end
        fprintf('No. of iterations needed:\t'); disp(jjj)
        
    case 2 %% l(T*z+t) = |Tz+t|_2 <= c'z+d
        disp('Lorentz cone')
        %% Initialiations  
        T = randn(3*n,2*n);   t = randn(3*n,1);
        c = randn(2*n,1); d = 3.9;
        K = [[T;c']; A];
        k = [[t;d]; b];        
        no.con.all = size(K,1);
        Mu = max(svds(K)); % Bound on K operator
        theta = 1; sigma = 1*.9/Mu; tau = (1/1)*sigma;
        p = zeros(3*n+1,1); old.p = p; tilde.p = p; % has the size of T + 1
        nu = zeros(n,1); old.nu = nu; tilde.nu = nu;
        z = zeros(2*n,1);   old.z = z;  bar.z = z;
        
        % precompute some stuff man
        matrix_inverse = inv(H+eye(size(H,1))*(1/tau));
        
        %% Main
        for jjj = 1:MAXITER
            %% step 1
            old.p = p;
            old.nu = nu;
            tilde.p = p + [sigma*(T*bar.z+t); sigma*(c'*z+d)];
            nu = nu + sigma*(A*bar.z-b);
           
            p = proj_secondOrderConeConj(tilde.p);

            %% step 2
            old.z = z;
            z = matrix_inverse * ((1/tau)*z-(K'*[p;nu]));

            %% step 3
            bar.z = z + theta*(z-old.z);

%             res.s  = tau*theta \ ( z-old.z);
%             res.r  = theta*K*(bar.z-old.z) + sigma \ ([p;nu]-[old.p;old.nu]);
            res.s  = 1/tau*theta*(z-old.z);
            res.r  = theta*K*(bar.z-old.z) + 1/sigma*([p;nu]-[old.p;old.nu]);
            
            rDual(jjj) = norm(res.s);    rPrimal(jjj) = norm(res.r);
            if (rDual(jjj) <= 1e-4 && rPrimal(jjj) <= 1e-4)
                    break;
            end
        end
        %%--------CVX---------%%
        cvx_begin
        cvx_solver 'sedumi'
            variables Z(2*n,1);
            obj = .5*Z'*H*Z;
            minimize ( obj )
            subject to 
            A*Z == b;
            { T*Z + t,c'*Z + d } == lorentz(3*n);
        cvx_end
        fprintf('Chambolle solution\n'); disp(z)
        fprintf('CVX solution\n'); disp(Z)
        fprintf('Inequality constraints active:\t');
        if (abs(norm(T*Z+t)-c'*Z-d)<=1e-8) 
            disp('YES') 
        else disp('NO') 
        end
        fprintf('No. of iterations needed:\t'); disp(jjj)
        
    case 3
        disp('l1-ball') %% l(T*z+t) = \delta_1(T*z+t,l), or |T*z+t|_1 <= l
        %% Initialiations
        T = randn(2*n);   t = randn(2*n,1);   l = 7.5;
        K = [T; A];
        k = [t; b];
        no.con.all = size(K,1);
        Mu = max(svds(K)); % Bound on K operator
        theta = 1; sigma = 1*.9/Mu; tau = (1/1)*sigma;
        p = zeros(2*n,1); old.p = p; tilde.p = p; % has the size of the primal z
        nu = zeros(n,1); old.nu = nu; tilde.nu = nu;
        z = zeros(2*n,1);   old.z = z;  bar.z = z;
        
        % precompute some stuff man
        matrix_inverse = inv(H+eye(size(H,1))*(1/tau));
        
        %% Main
        for jjj = 1:MAXITER 
            %% step 1
            old.p = p;
            old.nu = nu;
            tilde.p = p + sigma*(T*bar.z+t);
            nu = nu + sigma*(A*bar.z-b);

            p = prox_norm(tilde.p,l*sigma,inf);

            %% step 2
            old.z = z;
            z = matrix_inverse * ((1/tau)*z-(K'*[p;nu]));

            %% step 3
            bar.z = z + theta*(z-old.z);

%             res.s  = tau*theta \ ( z-old.z);
%             res.r  = theta*K*(bar.z-old.z) + sigma \ ([p;nu]-[old.p;old.nu]);
            res.s  = 1/tau*theta*( z-old.z);
            res.r  = theta*K*(bar.z-old.z) + 1/sigma*([p;nu]-[old.p;old.nu]);
            
            rDual(jjj) = norm(res.s);    rPrimal(jjj) = norm(res.r);
            if (rDual(jjj) <= 1e-3 && rPrimal(jjj) <= 1e-3)
                    break;
            end
        end
        %%--------CVX---------%%
        cvx_begin
        cvx_solver 'sedumi'
            variables Z(2*n,1);
            obj = .5*Z'*H*Z;
            minimize ( obj )
            subject to 
            A*Z == b;
            norm(T*Z+t,1) <= l;
        cvx_end
        fprintf('Chambolle solution\n'); disp(z)
        fprintf('CVX solution\n'); disp(Z)
        fprintf('Inequality constraints active:\t');
        if (abs(norm(T*Z+t,1)-l)<=1e-8) 
             disp('YES') 
        else disp('NO') 
        end
        fprintf('No. of iterations needed:\t'); disp(jjj)
        
    case 4
        disp('l2-ball') %% l(T*z+t) = \delta_2(T*z+t,l), or |T*z+t|_2 <= l
        %% Initialiations
        T = randn(2*n);   t = randn(2*n,1);   l = 3.5;
        K = [T; A];
        k = [t; b];
        no.con.all = size(K,1);
        Mu = max(svds(K)); % Bound on K operator
        theta = 1; sigma = 1*.9/Mu; tau = (1/1)*sigma;
        p = zeros(2*n,1); old.p = p; tilde.p = p; % has the size of the primal z
        nu = zeros(n,1); old.nu = nu; tilde.nu = nu;
        z = zeros(2*n,1);   old.z = z;  bar.z = z;
        
        % precompute some stuff man
        matrix_inverse = inv(H+eye(size(H,1))*(1/tau));
        
        %% Main
        for jjj = 1:MAXITER 
            %% step 1
            old.p = p;
            old.nu = nu;
            tilde.p = p + sigma*(T*bar.z+t);
            nu = nu + sigma*(A*bar.z-b);

            p = prox_norm(tilde.p,l*sigma,2);

            %% step 2
            old.z = z;
            z = matrix_inverse * ((1/tau)*z-(K'*[p;nu]));

            %% step 3
            bar.z = z + theta*(z-old.z);

%             res.s  = tau*theta \ ( z-old.z);
%             res.r  = theta*K*(bar.z-old.z) + sigma \ ([p;nu]-[old.p;old.nu]);
            res.s  = 1/tau*theta*( z-old.z);
            res.r  = theta*K*(bar.z-old.z) + 1/sigma*([p;nu]-[old.p;old.nu]);
            rDual(jjj) = norm(res.s);    rPrimal(jjj) = norm(res.r);
            if (rDual(jjj) <= 1e-3 && rPrimal(jjj) <= 1e-3)
                    break;
            end
        end
        %%--------CVX---------%%
        cvx_begin
        cvx_solver 'sedumi'
            variables Z(2*n,1);
            obj = .5*Z'*H*Z;
            minimize ( obj )
            subject to 
            A*Z == b;
            norm(T*Z+t,2) <= l;
        cvx_end
        fprintf('Chambolle solution\n'); disp(z)
        fprintf('CVX solution\n'); disp(Z)
        fprintf('Inequality constraints active:\t');
        if (abs(norm(T*Z+t,2)-l)<=1e-8) 
             disp('YES') 
        else disp('NO') 
        end
        fprintf('No. of iterations needed:\t'); disp(jjj)
        
    case 5
        disp('linfty-ball') %% l(T*z+t) = \delta_inf(T*z+t,l), or |T*z+t|_inf <= l
        %% Initialiations
        T = randn(2*n);   t = randn(2*n,1);   l = 3.3;
        K = [T; A];
        k = [t; b];
        no.con.all = size(K,1);
        Mu = max(svds(K)); % Bound on K operator
        theta = 1; sigma = 1*.9/Mu; tau = (1/1)*sigma;
        p = zeros(2*n,1); old.p = p; tilde.p = p; % has the size of the primal z
        nu = zeros(n,1); old.nu = nu; tilde.nu = nu;
        z = zeros(2*n,1);   old.z = z;  bar.z = z;
        
        %% Main
        for jjj = 1:MAXITER 
            %% step 1
            old.p = p;
            old.nu = nu;
            tilde.p = p + sigma*(T*bar.z+t);
            nu = nu + sigma*(A*bar.z-b);

            p = prox_norm(tilde.p,l*sigma,1);

            %% step 2
            old.z = z;
            z = matrix_inverse * ((1/tau)*z-(K'*[p;nu]));

            %% step 3
            bar.z = z + theta*(z-old.z);

%             res.s  = tau*theta \ ( z-old.z);
%             res.r  = theta*K*(bar.z-old.z) + sigma \ ([p;nu]-[old.p;old.nu]);

            res.s  = 1/tau*theta*( z-old.z);
            res.r  = theta*K*(bar.z-old.z) + 1/sigma*([p;nu]-[old.p;old.nu]);
            rDual(jjj) = norm(res.s);    rPrimal(jjj) = norm(res.r);
            if (rDual(jjj) <= 1e-3 && rPrimal(jjj) <= 1e-3)
                    break;
            end
        end
        %%--------CVX---------%%
        cvx_begin
        cvx_solver 'sedumi'
            variables Z(2*n,1);
            obj = .5*Z'*H*Z;
            minimize ( obj )
            subject to 
            A*Z == b;
            norm(T*Z+t,inf) <= l;
        cvx_end
        fprintf('Chambolle solution\n'); disp(z)
        fprintf('CVX solution\n'); disp(Z)
        fprintf('Inequality constraints active:\t');
        if (abs(norm(T*Z+t,inf)-l)<=1e-8) 
             disp('YES') 
        else disp('NO') 
        end
        fprintf('No. of iterations needed:\t'); disp(jjj)
        
    case 6
        disp('l1-norm') %% l*l(T*z+t) = l*|T*z+t|_1
        %% Initialiations
        T = randn(2*n);   t = randn(2*n,1);   l = 3.3;
        K = [T; A];
        k = [t; b];
        no.con.all = size(K,1);
        Mu = max(svds(K)); % Bound on K operator
        theta = 1; sigma = 1*.9/Mu; tau = (1/1)*sigma;
        p = zeros(2*n,1); old.p = p; tilde.p = p; % has the size of the primal z
        nu = zeros(n,1); old.nu = nu; tilde.nu = nu;
        z = zeros(2*n,1);   old.z = z;  bar.z = z;
        
        %% Main
        for jjj = 1:MAXITER 
            %% step 1
            old.p = p;
            old.nu = nu;
            tilde.p = p + sigma*(T*bar.z+t);
            nu = nu + sigma*(A*bar.z-b);

            % p = proj_box(tilde.p, -l*ones(2*n,1), l*ones(2*n,1));
            p = proj_normBall(tilde.p, Inf, l);

            %% step 2
            old.z = z;
            z = matrix_inverse * ((1/tau)*z-(K'*[p;nu]));

            %% step 3
            bar.z = z + theta*(z-old.z);

%             res.s  = tau*theta \ ( z-old.z);
%             res.r  = theta*K*(bar.z-old.z) + sigma \ ([p;nu]-[old.p;old.nu]);
            res.s  = 1/tau*theta*( z-old.z);
            res.r  = theta*K*(bar.z-old.z) + 1/sigma*([p;nu]-[old.p;old.nu]);
            
            rDual(jjj) = norm(res.s);    rPrimal(jjj) = norm(res.r);
            if (rDual(jjj) <= 1e-3 && rPrimal(jjj) <= 1e-3)
                    break;
            end
        end
        %%--------CVX---------%%
        cvx_begin
        cvx_solver 'sedumi'
            variables Z(2*n,1);
            obj = .5*Z'*H*Z + l*norm(T*Z+t,1);
            minimize ( obj )
            subject to 
            A*Z == b;
        cvx_end
        fprintf('Chambolle solution\n'); disp(z)
        fprintf('CVX solution\n'); disp(Z)
        fprintf('No. of iterations needed:\t'); disp(jjj)
        
    case 7
        disp('l2-norm') %% l*l(T*z+t) = l*|T*z+t|_2
        %% Initialiations
        T = randn(2*n);   t = randn(2*n,1);   l = 5.5;
        K = [T; A];
        k = [t; b];
        no.con.all = size(K,1);
        Mu = max(svds(K)); % Bound on K operator
        theta = 1; sigma = 1*.9/Mu; tau = (1/1)*sigma;
        p = zeros(2*n,1); old.p = p; tilde.p = p; % has the size of the primal z
        nu = zeros(n,1); old.nu = nu; tilde.nu = nu;
        z = zeros(2*n,1);   old.z = z;  bar.z = z;
        
        %% Main
        for jjj = 1:MAXITER 
            %% step 1
            old.p = p;
            old.nu = nu;
            tilde.p = p + sigma*(T*bar.z+t);
            nu = nu + sigma*(A*bar.z-b);

            p = proj_normBall(tilde.p, 2, l);

            %% step 2
            old.z = z;
            z = matrix_invers * ((1/tau)*z-(K'*[p;nu]));

            %% step 3
            bar.z = z + theta*(z-old.z);

%             res.s  = tau*theta \ ( z-old.z);
%             res.r  = theta*K*(bar.z-old.z) + sigma \ ([p;nu]-[old.p;old.nu]);
            res.s  = 1/tau*theta*( z-old.z);
            res.r  = theta*K*(bar.z-old.z) + 1/sigma*([p;nu]-[old.p;old.nu]);
            
            rDual(jjj) = norm(res.s);    rPrimal(jjj) = norm(res.r);
            if (rDual(jjj) <= 1e-3 && rPrimal(jjj) <= 1e-3)
                    break;
            end
        end
        %%--------CVX---------%%
        cvx_begin
        cvx_solver 'sedumi'
            variables Z(2*n,1);
            obj = .5*Z'*H*Z + l*norm(T*Z+t,2);
            minimize ( obj )
            subject to 
            A*Z == b;
        cvx_end
        fprintf('Chambolle solution\n'); disp(z)
        fprintf('CVX solution\n'); disp(Z)
        fprintf('No. of iterations needed:\t'); disp(jjj)
        
    case 8
        disp('linfty-norm') %% l*l(T*z+t) = l*|T*z+t|_inf
        %% Initialiations
        T = randn(2*n);   t = randn(2*n,1);   l = 5.5;
        K = [T; A];
        k = [t; b];
        no.con.all = size(K,1);
        Mu = max(svds(K)); % Bound on K operator
        theta = 1; sigma = 1*.9/Mu; tau = (1/1)*sigma;
        p = zeros(2*n,1); old.p = p; tilde.p = p; % has the size of the primal z
        nu = zeros(n,1); old.nu = nu; tilde.nu = nu;
        z = zeros(2*n,1);   old.z = z;  bar.z = z;
        
        %% Main
        for jjj = 1:MAXITER 
            %% step 1
            old.p = p;
            old.nu = nu;
            tilde.p = p + sigma*(T*bar.z+t);
            nu = nu + sigma*(A*bar.z-b);

            p = proj_normBall(tilde.p, 1, l);

            %% step 2
            old.z = z;
            z = matrix_invers * ((1/tau)*z-(K'*[p;nu]));

            %% step 3
            bar.z = z + theta*(z-old.z);

%             res.s  = tau*theta \ ( z-old.z);
%             res.r  = theta*K*(bar.z-old.z) + sigma \ ([p;nu]-[old.p;old.nu]);
            res.s  = 1/tau*theta*( z-old.z);
            res.r  = theta*K*(bar.z-old.z) + 1/sigma*([p;nu]-[old.p;old.nu]);
            rDual(jjj) = norm(res.s);    rPrimal(jjj) = norm(res.r);
            if (rDual(jjj) <= 1e-3 && rPrimal(jjj) <= 1e-3)
                    break;
            end
        end
        %%--------CVX---------%%
        cvx_begin
        cvx_solver 'sedumi'
            variables Z(2*n,1);
            obj = .5*Z'*H*Z + l*norm(T*Z+t,inf);
            minimize ( obj )
            subject to 
            A*Z == b;
        cvx_end
        fprintf('Chambolle solution\n'); disp(z)
        fprintf('CVX solution\n'); disp(Z)
        fprintf('No. of iterations needed:\t'); disp(jjj)
        
end
