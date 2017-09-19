/* ========================================================================== */
/* === ldl.h:  include file for the LDL package ============================= */
/* ========================================================================== */

/* Copyright (c) Timothy A Davis, http://www.suitesparse.com.
 * All Rights Reserved.  See README for the License.
 */

#include "user_suiteSparse_config.h"
#include "user_foo_data.h"
#define real data_t_primal_out


#ifdef LDL_LONG
#define LDL_int SuiteSparse_long
#define LDL_ID SuiteSparse_long_id

#define LDL_symbolic ldl_l_symbolic
#define LDL_numeric ldl_l_numeric
#define LDL_lsolve ldl_l_lsolve
#define LDL_dsolve ldl_l_dsolve
#define LDL_ltsolve ldl_l_ltsolve
#define LDL_perm ldl_l_perm
#define LDL_permt ldl_l_permt
#define LDL_valid_perm ldl_l_valid_perm
#define LDL_valid_matrix ldl_l_valid_matrix

#define chol_lsolve chol_l_lsolve
#define chol_ltsolve chol_l_ltsolve

#else
#define LDL_int int
#define LDL_ID "%d"

#define LDL_symbolic ldl_symbolic
#define LDL_numeric ldl_numeric
#define LDL_lsolve ldl_lsolve
#define LDL_dsolve ldl_dsolve
#define LDL_ltsolve ldl_ltsolve
#define LDL_perm ldl_perm
#define LDL_permt ldl_permt
#define LDL_valid_perm ldl_valid_perm
#define LDL_valid_matrix ldl_valid_matrix

#define chol_lsolve chol_lsolve

#define chol_ltsolve chol_ltsolve

#endif

/* ========================================================================== */
/* === int version ========================================================== */
/* ========================================================================== */

void ldl_symbolic (int n, int Ap [ ], int Ai [ ], int Lp [ ],
    int Parent [ ], int Lnz [ ], int Flag [ ], int P [ ], int Pinv [ ]) ;

int ldl_numeric (int n, int Ap [ ], int Ai [ ], real Ax [ ],
    int Lp [ ], int Parent [ ], int Lnz [ ], int Li [ ], real Lx [ ],
    real D [ ], real Y [ ], int Pattern [ ], int Flag [ ],
    int P [ ], int Pinv [ ]) ;

void ldl_lsolve (int n, real X [ ], int Lp [ ], int Li [ ],
    real Lx [ ]) ;

void ldl_dsolve (int n, real X [ ], real D [ ]) ;

void ldl_ltsolve (int n, real X [ ], int Lp [ ], int Li [ ],
    real Lx [ ]) ;

/*void chol_lsolve (int n, real X [ ], int Lp [ ], int Li [ ],
    real Lx [ ]) ;

void chol_ltsolve (int n, real X [ ], int Lp [ ], int Li [ ],
    real Lx [ ]) ;
*/

void ldl_perm  (int n, real X [ ], real B [ ], int P [ ]) ;
void ldl_permt (int n, real X [ ], real B [ ], int P [ ]) ;

int ldl_valid_perm (int n, int P [ ], int Flag [ ]) ;
int ldl_valid_matrix ( int n, int Ap [ ], int Ai [ ]) ;

/* ========================================================================== */
/* === long version ========================================================= */
/* ========================================================================== */

void ldl_l_symbolic (SuiteSparse_long n, SuiteSparse_long Ap [ ],
    SuiteSparse_long Ai [ ], SuiteSparse_long Lp [ ],
    SuiteSparse_long Parent [ ], SuiteSparse_long Lnz [ ],
    SuiteSparse_long Flag [ ], SuiteSparse_long P [ ],
    SuiteSparse_long Pinv [ ]) ;

SuiteSparse_long ldl_l_numeric (SuiteSparse_long n, SuiteSparse_long Ap [ ],
    SuiteSparse_long Ai [ ], real Ax [ ], SuiteSparse_long Lp [ ],
    SuiteSparse_long Parent [ ], SuiteSparse_long Lnz [ ],
    SuiteSparse_long Li [ ], real Lx [ ], real D [ ], real Y [ ],
    SuiteSparse_long Pattern [ ], SuiteSparse_long Flag [ ],
    SuiteSparse_long P [ ], SuiteSparse_long Pinv [ ]) ;

void ldl_l_lsolve (SuiteSparse_long n, real X [ ], SuiteSparse_long Lp [ ],
    SuiteSparse_long Li [ ], real Lx [ ]) ;

void ldl_l_dsolve (SuiteSparse_long n, real X [ ], real D [ ]) ;

void ldl_l_ltsolve (SuiteSparse_long n, real X [ ], SuiteSparse_long Lp [ ],
    SuiteSparse_long Li [ ], real Lx [ ]) ;

/*void chol_l_lsolve (SuiteSparse_long n, real X [ ], SuiteSparse_long Lp [ ],
    SuiteSparse_long Li [ ], real Lx [ ]) ;

void chol_l_ltsolve (SuiteSparse_long n, real X [ ], SuiteSparse_long Lp [ ],
    SuiteSparse_long Li [ ], real Lx [ ]) ;
*/

void ldl_l_perm  (SuiteSparse_long n, real X [ ], real B [ ],
    SuiteSparse_long P [ ]) ;
void ldl_l_permt (SuiteSparse_long n, real X [ ], real B [ ],
    SuiteSparse_long P [ ]) ;

SuiteSparse_long ldl_l_valid_perm (SuiteSparse_long n, SuiteSparse_long P [ ],
    SuiteSparse_long Flag [ ]) ;
SuiteSparse_long ldl_l_valid_matrix ( SuiteSparse_long n,
    SuiteSparse_long Ap [ ], SuiteSparse_long Ai [ ]) ;

/* ========================================================================== */
/* === LDL version ========================================================== */
/* ========================================================================== */

#define LDL_DATE "Oct 10, 2014"
#define LDL_VERSION_CODE(main,sub) ((main) * 1000 + (sub))
#define LDL_MAIN_VERSION 2
#define LDL_SUB_VERSION 2
#define LDL_SUBSUB_VERSION 1
#define LDL_VERSION LDL_VERSION_CODE(LDL_MAIN_VERSION,LDL_SUB_VERSION)

