// Simple matrix operation library

// NOTE : Add return vectors are assumed to be allocated by caller
// NOTE : No memory violation checks are made

#ifndef __matrix_ops__
#define __matrix_ops__
#include <math.h>
#include <string.h>

#include "user_foo_data.h"
#define real data_t_primal_out


void print_vector(const char *name, const real *x, const int len);

// y = x
//#define copy_vector(y,x,len) memcpy(y,x,sizeof(real)*len)
void copy_vector( real *y, const real *x, int len);

// returns sum(x)
real sum_vector(const real *x, const int len);

// returns sum(x.*x)
real sum_squared(const real *x, const int len);

// returns sum(|x|)
real sum_absolute(const real *x, const int len);

// y = alpha*x
void scale_vector(real *y, const real alpha, const real *x, const int len);

// returns <x,y>
real dot_vector(const real *x, const real *y, const int len);
void vector_constant_addition(real *y, const real *x, real constant, const int len);
#if 0
void divide_vector_elements(const real *y, const real *x, const int len );

#endif
// y = M*x
void matvec_product(const real *M, const real *x, real *y, const int nrows, const int ncols);

// Macro to extend an operation elementwise
#define forall(range_h) for(i=0; i<range_h; i++)

// Proximal functions
void prox_norm_one(real *xprox, const real *x, const real t, const int len);
void prox_norm_two(real *xprox, const real *x, const real t, const int len);
//void prox_norm_inf(real *x, real *z, real rho);
//void proj_normball_one(real *x, const real *z, const real c);

void proj_normball_one(real *xproj, real *x, const real c, const int len);
// Projections
void proj_normball_two(real *xproj, const real *x, const real c, const int len);
void proj_normball_inf(real *xproj, const real *x, const real c, const int len);
void proj_secondOrderCone(real *xproj, const real *x, const int len); // Project the point x = [t;y] onto the second order cone
void proj_secondOrderCone_conj(real *xproj, const real *x, const int len);
void proj_negative(real *xproj, const real *x, const int len);
void proj_positive(real *xproj, const real *x, const int len);

// Not used
void proj_box(real *xproj, const real *x, const real *lb, const real *ub, const int len);

// Vector norms
real norm_one(const real *x, const int len);
real norm_two(const real *x, const int len);
real norm_inf(const real *x, const int len);

void quick_sort (real *a, int n);
#endif
