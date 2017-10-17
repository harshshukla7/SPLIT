//
//  main.h
//  chol_mod_harsh
//
//  Created by harsh shukla on 23/01/15.
//  Copyright (c) 2015 harsh shukla. All rights reserved.
//

#ifndef chol_mod_harsh_main_h
#define chol_mod_harsh_main_h


#endif

#include <stdio.h>

#include "Accelerate/Accelerate.h"

//////cholmod

#include "cholmod.h"

////////////ldl

#include "ldl.h"

#define N 10	/* A is 10-by-10 */
#define ANZ 19	/* # of nonzeros on diagonal and upper triangular part of A */
#define LNZ 13	/* # of nonzeros below the diagonal of L */




//int ldl_suitesparse (int Ap[],int Ai[],double Ax[],double b[],int Parent,int Flag, int Lp[], int Li[],double Lx[], double Lnz[], double D[], double Y[], int Pattern, int d, int i);

int ldl_suitesparse_prefactor (int Ap[],int Ai[],double Ax[],double D[], double Lx[], int Li[], int Lp[], int d[1]);

int ldl_suitesparse_solve(double b[], double D[], double Lx[], int Li[], int Lp[], int d[1]);

void printmatrix ( unsigned m , unsigned n, double A[m*n]);
void printvector (unsigned m, double x[m]);