//
//  ldl_suitesparse.c
//  chol_mod_harsh
//
//  Created by harsh shukla on 23/01/15.
//  Copyright (c) 2015 harsh shukla. All rights reserved.
//

#include "main.h"

int ldl_suitesparse_prefactor (int Ap[],int Ai[],double Ax[], double D[], double Lx[], int Li[], int Lp[], int d[1]){
    
    double Y[N];
    // double Y [N] ;
    int  Parent [N], Lnz [N], Flag [N], Pattern [N] ;
    
    
    ldl_symbolic (N, Ap, Ai, Lp, Parent, Lnz, Flag, NULL, NULL) ;
    printf ("Nonzeros in L, excluding diagonal: %d\n", Lp [N]) ;
    d[1] = ldl_numeric (N, Ap, Ai, Ax, Lp, Parent, Lnz, Li, Lx, D, Y, Pattern,
                     Flag, NULL, NULL) ;
    
    return d[1];
}

int ldl_suitesparse_solve(double b[], double D[], double Lx[], int Li[], int Lp[], int d[1]){

    if (d[1] == N)
    {
        /* solve Ax=b, overwriting b with the solution x */
        ldl_lsolve (N, b, Lp, Li, Lx) ;
        ldl_dsolve (N, b, D) ;
        ldl_ltsolve (N, b, Lp, Li, Lx) ;
        for (int i = 0 ; i < N ; i++) printf ("x [%d] = %g\n", i, b [i]) ;
    }
    else
    {
        printf ("ldl_numeric failed, D (%d,%d) is zero\n", d[1], d[1]) ;
    }
    
    
    //tp_return = ldl_suitesparse ( Ap, Ai, Ax, b, Parent, Flag, Lp, Li, Lx, Lnz, D, Y, Pattern);
    
    return 0;
}


void printmatrix ( unsigned m , unsigned n, double A[m*n]){
    
    
    for (unsigned i=0; i<m; i++) {
        
        printf("\n");
        
        for (unsigned j=0; j<n; j++) {
            
            printf("%f \t \t", A[(i*n)+j]);
        }
    }
    
    printf("\n");
}

void printvector (unsigned m, double x[m]){
    
    printf("\n");
    
    for (unsigned i=0; i<m; i++){
        
        printf("%f \t \t", x[i]);
        
    }
    printf("\n");
    
}

