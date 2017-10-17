#ifndef __splitData_h__
#define __splitData_h__
<<<<<<< HEAD
#include <stdio.h>
#include <string.h>
double x[10];
#define x_len 10
double sol[20];
#define sol_len 20
double Arow[30];
#define Arow_len 30
double Acol[30];
#define Acol_len 30
int I[30];
#define I_len 30
=======
double x[10];
#define x_len = 10
double sol[20];
#define sol_len = 20
double Arow[30];
#define Arow_len = 30
double Acol[30];
#define Acol_len = 30
int I[30];
#define I_len = 30


void loadData() {
 Data *dat = splitLoad("splitData.dat");
 copyVar(x, dat, "x");
 copyVar(sol, dat, "sol");
 copyVar(Arow, dat, "Arow");
 copyVar(Acol, dat, "Acol");
 copyVar(I, dat, "I");
freeData(dat);
}


#endif
>>>>>>> origin/feature/lapack_suitesparse
