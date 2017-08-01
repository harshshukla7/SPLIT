


method = 'chol'; %% chol or ldl

hsd=splitData;


%%%% Next we generate the problem
%%%% methods avaialable are following
%%%% 'ldl_lp' ldl using clapack
%%%% 'ldl_ss' ldl using suitesparse
%%%% 'ldl'    ldl using matlab (I have not checked it if it works with this m function)
%%%% 'invert' preinversion and saving the matrix (I have not checked it if it works with this m function)
switch method
    case 'ldl'
        
        A=[1, 2, 0 ; 2 3 0 ; 0 0 5];
        b=[11;17;10];
        a_name=sprintf('rhs');
        hsd.add_var(a_name, b, 'type', 'real');
        hsd.add_function(coderFunc_mldivide('custom_KKT_solve',A, 'method', 'ldl_ss'));
        
        LNZ=nnz(A)-nnz(triu(A));
        ANZ=nnz(triu(A));
        N=size(A,2);
        %%%% Do not delete the following ....
        
        hsd.define('nn_lp',length(b),'int');
        hsd.define('LNZ_ss', LNZ, 'int');
        hsd.define('ANZ_ss', ANZ, 'int');
        hsd.define('N_ss', N, 'int');
        
        hsd.hl('#define ldl\n');
        
    case 'chol'
        
        
        hsd.hl('#define chol\n');
        Lin_Solve = 'llt';
        
        %%% Random PD matrix
        
        A = randn(3);
        [U,ignore] = eig((A+A')/2); % (Don't really have to divide by 2)
        M = U*diag(abs(randn(3,1)))*U';
        
        %%%%
        
        b = [1; 2; 3];
        a_name=sprintf('rhs');
        hsd.add_var(a_name, b, 'type', 'real');
        chol_solve = coderFunc_chol_solve('custom_chol_solve', A*A', 'Astr', 'AA','method', Lin_Solve);
        hsd.add_function(chol_solve);
        
        if strcmp(chol_solve.desc,'ldl_ss')
        hsd.cl('double *Lx_ss;');
        hsd.cl('int *Li_ss;');
        hsd.hl('#define suitesparse_linsolve\n');
        hsd.hl('extern double *Lx_ss;');
        hsd.hl('extern int *Li_ss;');
        end
        
        
        K = A*A';
        LNZ=nnz(K)-nnz(triu(K));
        ANZ=nnz(triu(K));
        N=size(K,2);
        hsd.define('nn_lp',size(K,2),'int');
        hsd.define('LNZ_ss', LNZ, 'int');
        hsd.define('ANZ_ss', ANZ, 'int');
        hsd.define('N_ss', N, 'int');
        
        
end




%% Write the probData

hsd.write_to_file('probData');

%%

AA = sparse(A*A');

[L,s,p] = chol(AA,'lower','vector');

normest(L*L'-AA);
tp1 = L\b
tp2 = L'\tp1


system(' gcc -Wall -o lin_sol test_linsolv.c matrix_ops.c probData.c ldl.c splitLoad.c');
system('./lin_sol');

