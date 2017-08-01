#ifndef _PROJECTIONS__H__
#define _PROJECTIONS__H__
    
    // fix it once the box constraint could be recognised
    void proj_box(double *x, double *z, double *lb, double *ub, int size);

    void proj_negative(double *x, double *z, int size, double weight, char *norm_type, double constant);
    void proj_normBall(double *x, double *z, int size, double weight, char *norm_type, double constant);
    void proj_positive(double *x, double *z, int size, double weight, char *norm_type, double constant);
    void proj_quadball(double *x, double *z, int size, double weight, char *norm_type, double constant);
    void proj_socp(double *x, double *z, int size, double weight, char *norm_type, double constant);
    
    /* actually it is not a projection - but a computation of a prox
       operator of a norm - rename it later on*/ 
    void prox_norm(double *x, double *z, int size, double weight, char *norm_type, double constant);

    // I think I should add a separate subroutine for projection soc - conjugate

#endif