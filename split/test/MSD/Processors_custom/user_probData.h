/***********************************************
 * Header file automatically generated by SPLIT
 *
 * More information : la.epfl.ch
 ***********************************************/
#ifndef __user_probData_h__
#define __user_probData_h__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "foo_data.h"
#include "user_matrix_ops.h"
#define real data_t_primal_out
#define REAL data_t_primal_out
#define nParam 4
#define nPrimal 46
#define nDual 28
#define nEqCon 32
#define rho_init 0.050000
#define rho_inv_init 20.000000
#define prox_var_init 0.000000
extern real rho;

extern real rhoinv; 

extern real prox_var; 

#define nn_lp 78
#define LNZ_ss 200
#define ANZ_ss 246
#define N_ss 78

#define par_ex_rows 4
#define par_ex_cols 1
#define par_ex_len 4
extern real par_ex[4];

#define sol_x_rows 46
#define sol_x_cols 1
#define sol_x_len 46
extern real sol_x[46];

#define sol_lam_rows 28
#define sol_lam_cols 1
#define sol_lam_len 28
extern real sol_lam[28];

#define sol_y_rows 28
#define sol_y_cols 1
#define sol_y_len 28
extern real sol_y[28];

#define l_rows 28
#define l_cols 1
#define l_len 28
extern real l[28];

#define f_rows 46
#define f_cols 1
#define f_len 46
extern real f[46];

#define b_rows 32
#define b_cols 1
#define b_len 32
extern real b[32];

#define _var_1__b_rows 32
#define _var_1__b_cols 1
#define _var_1__b_len 32
extern real _var_1__b[32];

#define Ltrans_p_ss_rows 29
#define Ltrans_p_ss_cols 1
#define Ltrans_p_ss_len 29
extern int Ltrans_p_ss[29];

#define Ltrans_i_ss_rows 28
#define Ltrans_i_ss_cols 1
#define Ltrans_i_ss_len 28
extern int Ltrans_i_ss[28];

#define Ltrans_x_ss_rows 28
#define Ltrans_x_ss_cols 1
#define Ltrans_x_ss_len 28
extern real Ltrans_x_ss[28];

#define n_Ltrans_ss_rows 1
#define n_Ltrans_ss_cols 1
#define n_Ltrans_ss_len 1
extern int n_Ltrans_ss[1];

#define KKT_L_iDiag_rows 78
#define KKT_L_iDiag_cols 1
#define KKT_L_iDiag_len 78
extern real KKT_L_iDiag[78];

#define KKT_L_lower_tri_rows 78
#define KKT_L_lower_tri_cols 78
#define KKT_L_lower_tri_len 6084

#define KKT_L_lower_tri_nz_per_row_rows 1
#define KKT_L_lower_tri_nz_per_row_cols 78
#define KKT_L_lower_tri_nz_per_row_len 78
extern int KKT_L_lower_tri_nz_per_row[78];

#define KKT_L_lower_tri_ind_rows 1
#define KKT_L_lower_tri_ind_cols 340
#define KKT_L_lower_tri_ind_len 340
extern int KKT_L_lower_tri_ind[340];

#define KKT_L_lower_tri_dat_rows 1
#define KKT_L_lower_tri_dat_cols 340
#define KKT_L_lower_tri_dat_len 340
extern real KKT_L_lower_tri_dat[340];

#define KKT_L_p_rows 1
#define KKT_L_p_cols 78
#define KKT_L_p_len 78
extern int KKT_L_p[78];

#define KKT_L_s_rows 78
#define KKT_L_s_cols 1
#define KKT_L_s_len 78
extern real KKT_L_s[78];

#define _var_4__p_ss_rows 79
#define _var_4__p_ss_cols 1
#define _var_4__p_ss_len 79
extern int _var_4__p_ss[79];

#define _var_4__i_ss_rows 79
#define _var_4__i_ss_cols 1
#define _var_4__i_ss_len 79
extern int _var_4__i_ss[79];

#define _var_4__x_ss_rows 79
#define _var_4__x_ss_cols 1
#define _var_4__x_ss_len 79
extern real _var_4__x_ss[79];

#define n__var_4__ss_rows 1
#define n__var_4__ss_cols 1
#define n__var_4__ss_len 1
extern int n__var_4__ss[1];

#define KKT_LT_iDiag_rows 78
#define KKT_LT_iDiag_cols 1
#define KKT_LT_iDiag_len 78
extern real KKT_LT_iDiag[78];

#define KKT_LT_upper_tri_rows 78
#define KKT_LT_upper_tri_cols 78
#define KKT_LT_upper_tri_len 6084

#define KKT_LT_upper_tri_nz_per_row_rows 1
#define KKT_LT_upper_tri_nz_per_row_cols 78
#define KKT_LT_upper_tri_nz_per_row_len 78
extern int KKT_LT_upper_tri_nz_per_row[78];

#define KKT_LT_upper_tri_ind_rows 1
#define KKT_LT_upper_tri_ind_cols 340
#define KKT_LT_upper_tri_ind_len 340
extern int KKT_LT_upper_tri_ind[340];

#define KKT_LT_upper_tri_dat_rows 1
#define KKT_LT_upper_tri_dat_cols 340
#define KKT_LT_upper_tri_dat_len 340
extern real KKT_LT_upper_tri_dat[340];

#define KKT_LT_s_rows 78
#define KKT_LT_s_cols 1
#define KKT_LT_s_len 78
extern real KKT_LT_s[78];

#define L_p_ss_rows 47
#define L_p_ss_cols 1
#define L_p_ss_len 47
extern int L_p_ss[47];

#define L_i_ss_rows 28
#define L_i_ss_cols 1
#define L_i_ss_len 28
extern int L_i_ss[28];

#define L_x_ss_rows 28
#define L_x_ss_cols 1
#define L_x_ss_len 28
extern real L_x_ss[28];

#define n_L_ss_rows 1
#define n_L_ss_cols 1
#define n_L_ss_len 1
extern int n_L_ss[1];
#ifdef lapack_linsolve 
extern __CLPK_integer ipiv[nn_lp];
#endif 
extern double *Lx_ss;
extern int *Li_ss;


// Functions
void custom_compute_parametric(REAL l[28], REAL f[46], REAL b[32], REAL par[4]);
void compute_parametric_pB(REAL y[32], REAL x[4]);
void custom_mult_Ltrans(REAL y[46], REAL x[28]);
/*
ldl
*/
void custom_solve_kkt(REAL y[78], REAL x[78]);
void custom_solve_kkt_solve1(real _x[78], real _b[78]);
void custom_solve_kkt_solve2(REAL y[78], REAL x[78]);
void custom_solve_kkt_solve3(real _x[78], real _z[78]);
void custom_compute_prefactor();
void custom_mult_L(REAL workDual[28], REAL x[46]);
void custom_prox(real y[nDual], real x[nDual]);


#endif