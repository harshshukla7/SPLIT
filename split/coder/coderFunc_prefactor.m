classdef coderFunc_prefactor < coderFunc
    
   methods
       function  f=coderFunc_prefactor(func_name,A,method)
       
           %%% step 1: check the input arguments and switch to lapack or
           %%% suitesparse
           %%% step 2: select the proper name for the function           
           
           
           
           p = inputParser;
           addRequired(p, 'func_name',  @ischar);
           addRequired(p, 'A',          @isnumeric);
           addRequired(p,'method',      @ischar);           
           
           
           
           f = f@coderFunc('void %s()', func_name);
      
           switch method
               
               case 'lapack'
                   
                   %%%step 1 Read the matrix in lapack format : just vec
                   %%%operation :)
                   %%%step 2 write the variables 
                   %%%step 3 call the factor function
                   
                   
                   %%%% Step 1
                   
                   nn=max(size(A,1),size(A,2));
                   KKT_lp=vec(full(A'));
                   
                %   KKT_lp
                   %%%%% Step 2
                   
                   a_length=length(KKT_lp);
                   a_name=sprintf('KKT_lp');
                   f.add_var(a_name,KKT_lp,'type','real');
                   %f.define('nn_lp',nn, 'type', 'int');
                   
                   %%%% Step 3 
                   
                   
                   f.pl('char uplo=''U'';');
                   f.pl('__CLPK_integer info, n, lda ;');
                   f.pl('n=nn_lp;');
                   f.pl('lda=nn_lp;');
                   f.pl('__CLPK_integer lwork=-1;');
                   f.pl('__CLPK_doublereal *work=malloc(sizeof(__CLPK_doublereal));');
                   f.pl('dsytrf_(&uplo, &n, KKT_lp, &lda, ipiv, work, &lwork, &info);');
                   f.pl('lwork=work[0];');
                   f.pl('free(work);');
                   f.pl('work=malloc(lwork*sizeof(__CLPK_doublereal));');
                   f.pl('dsytrf_(&uplo, &n, KKT_lp, &lda, ipiv, work, &lwork, &info);');
                   f.pl('free(work);');
                   
                   %f.pl('ldl_lapack_prefactor(KKT_lp);');                  
                   
                   
                   
               case 'suitesparse'
                   
                   %%%step 1 Read  the matrix in suitesparse
                   %%%format so the function ;)
                   %%%Step 2: write the variable 
                   %%%step 3 write the factor function for
                   %%%suitesparse: fprintf ;)
                   
                   %%%% step 1
                   [Ap_ss,Ai_ss,Ax_ss] =suitesparse_format(A);
                   
                    LNZ=nnz(A)-nnz(triu(A));
                    ANZ=nnz(triu(A));
                    N=size(A,2);
                   %%%% step 2
                   
                   a_name=sprintf('Ap_ss');
                   f.add_var(a_name,Ap_ss, 'type', 'int');
                   %Ap_ss
                   
                   a_name=sprintf('Ai_ss');
                   f.add_var(a_name,Ai_ss, 'type', 'int');
                   %Ai_ss
                   
                   a_name=sprintf('Ax_ss');
                   f.add_var(a_name,Ax_ss, 'type', 'real');
                   %Ax_ss
                   
                   %a_name=sprintf('Lx_ss');                   
                   %f.add_var(a_name,zeros(LNZ,1),'type', 'real');
                   
                   
                   a_name=sprintf('D_ss');                   
                   f.add_var(a_name,zeros(N,1),'type', 'real');
                   
                   
                   %a_name=sprintf('Li_ss');                   
                   %f.add_var(a_name,zeros(LNZ,1),'type', 'int');
                   
                   a_name=sprintf('Lp_ss');                   
                   f.add_var(a_name,zeros((N+1),1),'type', 'int');
                   
                   a_name=sprintf('d_ss');                   
                   f.add_var(a_name,zeros(1,1),'type', 'int');
                   
                   
                   %%%% step 3
                   
                   f.pl('real Y[N_ss];');
                   f.pl('int  Parent [N_ss], Lnz [N_ss], Flag [N_ss], Pattern [N_ss];');
                   f.pl('ldl_symbolic (N_ss, Ap_ss, Ai_ss, Lp_ss, Parent, Lnz, Flag, NULL, NULL);');
                   f.pl('int tmp_ss=Lp_ss[N_ss];');
                   f.pl('Li_ss = malloc(sizeof(int)*tmp_ss);');
                   f.pl('Lx_ss = malloc(sizeof(real)*tmp_ss);');
                   %f.pl('printf ("Nonzeros in L, excluding diagonal: %%d\n", Lp [N_ss]) ;');
                   f.pl('d_ss[0] = ldl_numeric (N_ss, Ap_ss, Ai_ss, Ax_ss, Lp_ss, Parent, Lnz, Li_ss, Lx_ss, D_ss, Y, Pattern, Flag, NULL, NULL) ;');
                   %f.pl('printf ("d_ss is : %%d\n", d_ss[0]) ;');
                   
                   %f.pl('ldl_suitesparse_prefactor ( Ap_ss, Ai_ss, Ax_ss, D_ss, Lx_ss, Li_ss, Lp_ss, d_ss);');
                   
                   
                   
                                      
               case 'none'
                   
                   f.pl('printf("No prefactorization in c is needed for this method");');
           end
       end
   end
end