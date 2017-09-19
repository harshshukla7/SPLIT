//
//  main.c
//  chol_mod_harsh
//
//  Created by harsh shukla on 05/01/15.
//  Copyright (c) 2015 harsh shukla. All rights reserved.
//

#include "main.h"

int main(int argc, const char * argv[]) {
    
    
    
    ////////////ldl Demostration starts here
    
        //int tp =10;
    //double *Y;
    //Y = (double*)malloc(tp*sizeof(double));
    
    //double *tp2;
    //memset(tp2, 0, tp*sizeof(double));
    
    int    Ap [N+1] = {0, 1, 2, 3, 4,   6, 7,   9,   11,      15,     ANZ},
    Ai [ANZ] = {0, 1, 2, 3, 1,4, 5, 4,6, 4,7, 0,4,7,8, 1,4,6,9 } ;
    double Ax [ANZ] = {1.7, 1., 1.5, 1.1, .02,2.6, 1.2, .16,1.3, .09,1.6,
        .13,.52,.11,1.4, .01,.53,.56,3.1},
    b [N] = {.287, .22, .45, .44, 2.486, .72, 1.55, 1.424, 1.621, 3.759};
    double Lx [LNZ], D [N], Y[N];
    // double Y [N] ;
    int Li [LNZ], Lp [N+1], Parent [N], Lnz [N], Flag [N], Pattern [N], d[1] ;
    

    
    int tp_return;
    
    ////////////// start comment here
    
    /* factorize A into LDL' (P and Pinv not used) */
    /*
    ldl_symbolic (N, Ap, Ai, Lp, Parent, Lnz, Flag, NULL, NULL) ;
    printf ("Nonzeros in L, excluding diagonal: %d\n", Lp [N]) ;
    d = ldl_numeric (N, Ap, Ai, Ax, Lp, Parent, Lnz, Li, Lx, D, Y, Pattern,
                     Flag, NULL, NULL) ;
    
    if (d == N)
    {
                ldl_lsolve (N, b, Lp, Li, Lx) ;
        ldl_dsolve (N, b, D) ;
        ldl_ltsolve (N, b, Lp, Li, Lx) ;
        for (i = 0 ; i < N ; i++) printf ("x [%d] = %g\n", i, b [i]) ;
    }
    else
    {
        printf ("ldl_numeric failed, D (%d,%d) is zero\n", d, d) ;
    }
    
     ///////////// end comment here */
    
    tp_return = ldl_suitesparse_prefactor ( Ap, Ai, Ax, D, Lx, Li, Lp, d);
    tp_return = ldl_suitesparse_solve(b, D, Lx, Li, Lp, d);
    //////////////// ldl demostration ends here
    
    
    printvector(10, b);
    
    return (0) ;

    //static const N=10;
    //////////////// Cholmod demostration
    
    /*
     cholmod_sparse *A ;
     cholmod_dense *x, *b, *r ;
     cholmod_factor *L ;
     
     double one [2] = {1,0}, m1 [2] = {-1,0} ;
     cholmod_common c ;
     cholmod_start (&c) ;
     
     A = cholmod_read_sparse (stdin, &c) ;
     cholmod_print_sparse (A, "A", &c) ;
     if (A == NULL || A->stype == 0)
     {
     cholmod_free_sparse (&A, &c) ;
     cholmod_finish (&c) ;
     return (0) ;
     }
     
     
     b = cholmod_ones (A->nrow, 1, A->xtype, &c) ;
     L = cholmod_analyze (A, &c) ;
     cholmod_factorize (A, L, &c) ;
     x = cholmod_solve (CHOLMOD_A, L, b, &c) ;
     r = cholmod_copy_dense (b, &c) ;
     cholmod_sdmult (A, 0, m1, one, x, r, &c) ;
     printf ("norm(b-Ax) %8.1e\n",
     cholmod_norm_dense (r, 0, &c)) ;
     cholmod_free_factor (&L, &c) ;
     cholmod_free_sparse (&A, &c) ;
     cholmod_free_dense (&r, &c) ;
     cholmod_free_dense (&x, &c) ;
     cholmod_free_dense (&b, &c) ;
     cholmod_finish (&c) ;
     
     */
    
    ////////cholmod ends here

}
