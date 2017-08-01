classdef coderFunc_chol_solve < coderFunc
    
    methods
        function f = coderFunc_chol_solve(func_name, A, varargin)
            %
            % func_mldivide(dat, func_name, A, param1, val1, etc)
            %
            % Generate code to solve a system of equations
            %
            % y = A\x
            %
            % Function prototype:
            %   void func_name(y[n], x[n])
            %
            % Arguments:
            %  func_name  : string
            %  A          : Matrix data
            %
            % Optional:
            %  Astr       : Name to use when storing the matrix [otherwise random]
            %  method     : One of 'ldl', 'invert'
            %               If not specified, method will be chosen based on the
            %               structure of A
            %  desc       : Text description of the function
            %  zero_tol   : Value below which a number is considered zero [1e-10]
            %               Used to sparsify the invert method
            %
            
            p = inputParser;
            addRequired(p, 'func_name',  @ischar);
            addRequired(p, 'A',          @isnumeric);
            
            addParameter(p, 'Astr',       uniqueVarName, @ischar);
            addParameter(p, 'method', 'auto', ...
                @(x) any(validatestring(x,{'auto','invert','llt','ldl_ss','llt_ss','ldl_lp'})));
            addParameter(p, 'desc', '', @ischar);
            addParameter(p, 'zero_tol', 1e-10, @isnumeric);
            
            % Parse and copy the results to the workspace
            parse(p, func_name, A, varargin{:});
            cellfun(@(q) evalin('caller',[q ' = p.Results.' q ';']), fieldnames(p.Results))
            
            % Write the function prototype
            [m,n] = size(A);
            assert(m==n, 'func_mldivide currently only solves square matrices')
            f = f@coderFunc('void %s(REAL y[%i], REAL x[%i])', func_name, n, n);
            
            
            
            if strcmpi(method, 'auto')
                % Try and guess the right method
                warning('Harsh / Giorgos : Put code here to guess which approach should be used by analysing the matrix A and the timing plots')
                method = 'ldl_ss';
            end
            
            f.desc = method;
            
            % Write the code to actually solve the problem
            switch method
                %%
                case 'invert'
                    % Invert the matrix and solve the problem via matrix-vector
                    % multiplication
                    
                    
                    invA = inv(full(A));
                    invA(abs(invA) < zero_tol) = 0;
                    
                    
                    %                     figure()
                    %                     subplot(2,1,1)
                    %                     spy(invA)
                    %                     subplot(2,1,2)
                    %                     spy(A)
                    %
                    fInv = coderFunc_times([func_name '_inv'], invA, 'Astr', Astr,'method', 'auto');
                    
                    % Copy to the body of the function into ours
                    f.p(fInv.body);
                    
                    % Copy the data over
                    f.data = fInv.data;
                    
                    %%%%%%% Harsh edits here for including prefactor function which
                    %%%%%%% does nothing
                    
                    %f.add_func(coderFunc_prefactor('custom_compute_prefactor',A,'none'));
                    
                    
                    %%%%% Edited part from harsh ends here
                    
                    %%
                case 'llt'
                    
                    
                    
%                     if(~ispdmat(A*A'))
%                         disp('Matrix must be positive definite for cholesky factorization')
%                     end
                    
                    % Pre-factor the matrix A
                    A = sparse(A);
                    
                    [L,s,p] = chol(A,'lower','vector');
                    
                    if(s~=0)
                        error('Cholesky decomp failed')
                    end
                    
                    %                     [L,D,p,S] = ldl(A,0.01,'vector');
                    %                     s = diag(S);
                    %                     s = s(p);
                    
                    % Create all the functions that we will need:
                    
                    % y  = L \ (s .* b(p));
                    f.add_func(...
                        coderFunc_lower_triangular_solve_chol(L, [Astr '_L'], p, [func_name '_solve1']));
                    
                    %                     % z  = inv(D) * y;
                    %                     f.add_func(...
                    %                         coderFunc_times([func_name '_solve2'], inv(D), 'method', 'ss'));
                    %
                    % q  = L' \ z;
                    % x(p,1)  = s .* q;
                    f.add_func(...
                        coderFunc_upper_triangular_solve_chol(L', [Astr '_LT'], p, [func_name '_solve2']));
                    
                    % Create the function that ties everything together
                    f.pl('  %s_solve1(y,  x);', func_name)
                    
                    f.pl('  %s_solve2(x,  y);', func_name)
                    
                    
                    %%%%%%% Harsh edits here for including prefactor function which
                    %%%%%%% does nothing
                    
                    %f.add_func(coderFunc_prefactor('custom_compute_prefactor',A,'none'));
                    
                    %%%%% Edited part from harsh ends here
                    
                case 'llt_ss'
                    
%                     if(~ispdmat(A*A'))
%                         disp('Matrix must be positive definite for cholesky factorization')
%                     end
%                     
                    %f.add_func(coderFunc_prefactor('custom_compute_prefactor',A,'suitesparse'));
                    
                  
                    error('Not Implemented Yet');
                    A = sparse(A);
                    
                    [L,s,p] = chol(A,'lower','vector');
                    
                    if(s~=0)
                        error('Cholesky decomp failed')
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%% save L matrix 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    [Ap_ss,Ai_ss,Ax_ss] =sparse_to_csr(sparse(tril(L,-1)));
                    Ap_ss = Ap_ss-1;
                    Ai_ss = Ai_ss - 1;
                    Ai_ss = [0;Ai_ss];
                    %%% step 2
                   
                   a_name=sprintf('Lp_chol');
                   f.add_var(a_name,Ap_ss, 'type', 'int');
                   %Ap_ss
                   
                   a_name=sprintf('Li_chol');
                   f.add_var(a_name,Ai_ss, 'type', 'int');
                   %Ai_ss
                   
                   a_name=sprintf('Lx_chol');
                   f.add_var(a_name,Ax_ss, 'type', 'real');
                   %Ax_ss
                   
                   a_name=sprintf('p_perm');
                   f.add_var(a_name,full(p), 'type', 'int');
                   
                   a_name=sprintf('L_diag');
                   f.add_var(a_name,1./diag(L),'type','double');
                   %a_name=sprintf('Lx_ss');                   
                   %f.add_var(a_name,zeros(LNZ,1),'type', 'real');
                   
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    f.pl('copy_vector(y,x,%i);', n); %%% copy the KKTRHP
                    f.pl('chol_lsolve (%i, y, Lp_chol, Li_chol, Lx_chol) ;',n);
                    f.pl('chol_ltsolve (%i, y, Lp_chol, Li_chol, Lx_chol) ;',n);
                    %f.pl('copy_vector(x,y,%i);', n); %%% copy the KKTRHP
                    %f.pl('ldl_permt (%i, y, x, p_perm);',n);
                    
               case 'ldl_ss'
                    
%                     assert(norm(full(A - A')) < 1e-12, 'Matrix must be symmetric to use LDL factorization')
%                     
                    f.add_func(coderFunc_prefactor('custom_compute_prefactor',A,'suitesparse'));
                    
                    
                    f.pl('copy_vector(y,x,%i);', n); %%% copy the KKTRHP
                    %f.pl('ldl_suitesparse_solve(y, D_ss, Lx_ss, Li_ss, Lp_ss, d_ss);'); %%% call the solve function
                    
                    f.pl('if (d_ss[0] == N_ss){');
                    f.pl('ldl_lsolve (N_ss, y, Lp_ss, Li_ss, Lx_ss) ;');
                    f.pl('ldl_dsolve (N_ss, y, D_ss) ;');
                    f.pl('ldl_ltsolve (N_ss, y, Lp_ss, Li_ss, Lx_ss) ;');
                    %   f.pl('//char *sol ="solution_vector";');
                    %    f.pl('//int tp=N_ss;');
                    %     f.pl('//print_vector(sol, y, tp);');
                    f.pl('}');
                    %f.pl('for (int i = 0 ; i < N_ss ; i++) printf ("x [%d] = %g\n", i, x[i]) ;');
                    %f.pl('}');
                    f.pl('else');
                    f.pl('{');
                    f.pl(' printf ("------------- ATTENTION ---------------- ");');
                    f.pl('printf ("------------- ATTENTION ----------------");');
                    f.pl('printf ("ldl_numeric failed, D has digonal zeros.");');
                    f.pl('printf ("Factorization using SuiteSparse failed. We reccomend you to pre-inver the KKT matrix or do LDL factorization using MATLAB.");');
                    f.pl('char err_ss;');
                    f.pl('printf("Type a character to coninue.");');
                    f.pl('err_ss = getchar();');
                    %f.pl('}');
                    f.pl('}');
                     
                    
                    case 'ldl_lp'
                    
                    %%%% Step1 : check if matrix is symmetric
                    assert(norm(full(A - A')) < 1e-12, 'Matrix must be symmetric to use LDL factorization')
                    
                    
                    %%%% Step2 :  Add coder_func_prefactor for ....
                    %%% (i) adding prefector c function in header
                    %%% and c file but not in this function
                    %%% (ii) Edit this function to solve the lienar
                    %%% system
                    
                    
                    f.add_func(coderFunc_prefactor('custom_compute_prefactor',A,'lapack')); %%% adding prefactor function
                    
                    f.pl('copy_vector(y,x,%i);', n); %%% copy the KKTRHP
                    %f.pl('ldl_lapack_solve(KKT_lp, y);'); %%% call the solve function
                    
                    f.pl('__CLPK_integer nrhs=1;');
                    f.pl('__CLPK_integer ldb=nn_lp;');
                    f.pl('__CLPK_integer info, n, lda ;');
                    f.pl('n=nn_lp;');
                    f.pl('lda=nn_lp;');
                    f.pl('char uplo=''U'';');
                    f.pl('dsytrs_(&uplo, &n, &nrhs, KKT_lp, &lda, ipiv, y, &ldb, &info);');
                    %  f.pl('//char *sol ="solution_vector";');
                    %  f.pl('//int tp=nn_lp;');
                    %  f.pl('//print_vector(sol, y, tp);');
                
                    
            end
        end
    end
end


