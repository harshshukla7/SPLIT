classdef coderFunc_mldivide < coderFunc
    
    methods
        function f = coderFunc_mldivide(func_name, A, varargin)
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
            %  nPrimal    : Number of primal variables
            %  paral      : Level of parallelism (only for FPGAs, SoC)
            %  lat        : adder_latency only for FPGA
            
            
            p = inputParser;
            addRequired(p, 'func_name',  @ischar);
            addRequired(p, 'A',          @isnumeric);
            
            addParameter(p, 'Astr',       uniqueVarName, @ischar);
            addParameter(p, 'method', 'auto', ...
                @(x) any(validatestring(x,{'auto', 'auto_FPGA', 'invert','invert_FPGA', 'invert_FPGA_MAC', 'invert_FPGA_tree', 'ldl','ldl_lp','ldl_ss'})));
            addParameter(p, 'desc', '', @ischar);
            addParameter(p, 'zero_tol', 1e-10, @isnumeric);
            addParameter(p, 'nPrimal', @isnumeric);
            
            addParameter(p, 'lat', 8, @isnumeric);
            addParameter(p, 'paral',  @isnumeric);
            addParameter(p, 'hard', 2,  @isnumeric);
            
            % Parse and copy the results to the workspace
            parse(p, func_name, A, varargin{:});
            cellfun(@(q) evalin('caller',[q ' = p.Results.' q ';']), fieldnames(p.Results))
            
            % Write the function prototype
            [m,n] = size(A);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % TODO: guess the best level of parallelism for FPGA and SoC 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            assert(m==n, 'func_mldivide currently only solves square matrices')
            f = f@coderFunc('void %s(REAL y[%i], REAL x[%i])', func_name, n, n);
            
            
            
            if strcmpi(method, 'auto_FPGA')
                % Try and guess the right method
                
                method = 'invert_FPGA_MAC';
            end
            
            
             if strcmpi(method, 'auto')
                % Try and guess the right method
                
                method = 'ldl_ss';
             end
            
            f.desc = method;
            
            % Write the code to actually solve the problem
            switch method
                %%
                case 'invert_FPGA_MAC'
                    
                    %first invert the matrix
                    invA = inv(full(A));
%                     invA = invA(1:nPrimal,:);
                    
                    %do not save the matrix as it will be saved in c file
                    
                    % to do 
                    % 1. par_request parsing
                    % 2. adder latency in settings
                    
                    set_fpga.adder_lat = lat;
                    set_fpga.hard = hard;
                    parall = paral;
                    split_MV_MAC(invA, parall, set_fpga);
                    % call the mat-vec generator
                    
                    f.pl('mv_mult(y, x);', size(invA,1), size(invA,2) );
                    
                case 'invert_FPGA_tree'
                    
                    set_fpga.adder_lat = lat;
                    set_fpga.hard = hard;
                    split_MV_MAC(H, paral, set_fpga)
                    split_MV_tree(H, settings)
                    
                    % call the mat-vec generator
                    
                    f.pl('mv_mult(y, x);', size(invA,1), size(invA,2) );
                    
                    
                case 'invert_FPGA'
                    
                    invA = inv(full(A));
                    invA = invA(1:nPrimal,:);
                    
                    f.add_var(Astr, invA');
                    
                     f.pl('int i, j, i_offset;', n); 
                     f.pl('reset_loop_kkt: for(i = 0; i < %d; i++){', size(invA,1));
                     f.pl('  #pragma HLS PIPELINE');
                     f.pl('  y[i] = 0.0;');
                     f.pl('}');
                     f.pl('col_loop: for(i = 0;i < %d; i++){', size(invA,2));
                     f.pl('  row_loop: for(j=0 ,  i_offset = 0; j < %d; j++, i_offset += %d){', size(invA,1), size(invA,2));
                     f.pl('#pragma HLS LOOP_FLATTEN');
                     f.pl('#pragma HLS PIPELINE');
                     f.pl('    y[j]=KKT[i_offset+i]*x[i] + y[j];');
                     f.pl('  }');
                     f.pl('}');
                     
                     fileID = fopen('user_mv_mult.cpp','w');
                     fclose(fileID);
                     
                     fileID = fopen('user_mv_mult.h','w');
                     fclose(fileID);
                    % f.pl('copy_vector(y,x,%i);', n); %%% copy the KKTRHP
                    
                    f.add_func(coderFunc_prefactor('custom_compute_prefactor',A,'none'));
                    
                
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
                    fInv = coderFunc_times([func_name '_inv'], invA, 'Astr', Astr);
                    
                    % Copy to the body of the function into ours
                    f.p(fInv.body);
                    
                    % Copy the data over
                    f.data = fInv.data;
                    
                    %%%%%%% Harsh edits here for including prefactor function which
                    %%%%%%% does nothing
                    
                    f.add_func(coderFunc_prefactor('custom_compute_prefactor',A,'none'));
                    
                    
                    %%%%% Edited part from harsh ends here
                    
                    %%
                case 'ldl'
                    
                    
                    % Solve the symmetric system of equations via LDL pre-factorization
                    %assert(norm(full(A - A')) == 0, 'Matrix must be symmetric to use LDL factorization')
                    
                   assert(norm(full(A - A')) < 1e-12, 'Matrix must be symmetric to use LDL factorization')
                    % Pre-factor the matrix A
                    A = sparse(A);
                    [L,D,p,S] = ldl(A,0.01,'vector');
                    s = diag(S);
                    s = s(p);
                    
                    % Create all the functions that we will need:
                    
                    % y  = L \ (s .* b(p));
                    f.add_func(...
                        coderFunc_lower_triangular_solve(L, [Astr '_L'], s, p, [func_name '_solve1']));
                    
                    % z  = inv(D) * y;
                    f.add_func(...
                        coderFunc_times([func_name '_solve2'], inv(D), 'method', 'ss'));
                    
                    % q  = L' \ z;
                    % x(p,1)  = s .* q;
                    f.add_func(...
                        coderFunc_upper_triangular_solve(L', [Astr '_LT'], s, p, [func_name '_solve3']));
                    
                    % Create the function that ties everything together
                    f.pl('  %s_solve1(y,  x);', func_name)
                    f.pl('  %s_solve2(x,  y);', func_name)
                    f.pl('  %s_solve3(y,  x);', func_name)
                    
                    
                    %%%%%%% Harsh edits here for including prefactor function which
                    %%%%%%% does nothing
                    
                    f.add_func(coderFunc_prefactor('custom_compute_prefactor',A,'none'));
                    
                    %%%%% Edited part from harsh ends here
                    
                    
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
                    
                    
                case 'ldl_ss'
                    
                    assert(norm(full(A - A')) < 1e-12, 'Matrix must be symmetric to use LDL factorization')
                    
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
                    
                 case 'chol'
                     
                     f.add_func(coderFunc_prefactor('custom_compute_prefactor',A,'none'));
                     
                    % step 1: set vector entries to zeros
                    % step 2: set another vector entries to zeros
                    % step 3: mat-vec
                    % step 4: vector addition
                    % step 5: lower cholesky solve
                    % step 6: upper cholesky solve
                    % step 7: mat-vec
                    
                    
                    
            end
        end
    end
end


