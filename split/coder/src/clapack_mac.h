/*
   =================================================================================================
   Definitions and prototypes for LAPACK as provided Apple Computer.

   Documentation of the LAPACK interfaces, including reference implementations, can be found on the web
   starting from the LAPACK FAQ page at this URL (verified live as of April 2002):
        http://netlib.org/lapack/faq.html
        
   A hardcopy maanual is:
        LAPACK Users' Guide, Third Edition. 
        @BOOK{laug,
            AUTHOR = {Anderson, E. and Bai, Z. and Bischof, C. and
                        Blackford, S. and Demmel, J. and Dongarra, J. and
                        Du Croz, J. and Greenbaum, A. and Hammarling, S. and
                        McKenney, A. and Sorensen, D.},
            TITLE = {{LAPACK} Users' Guide},
            EDITION = {Third},
            PUBLISHER = {Society for Industrial and Applied Mathematics},
            YEAR = {1999},
            ADDRESS = {Philadelphia, PA},
            ISBN = {0-89871-447-8 (paperback)} }

   =================================================================================================
*/
#ifndef __CLAPACK_H
#define __CLAPACK_H
 
#ifdef __cplusplus
extern "C" {
#endif


#if defined(__LP64__) /* In LP64 match sizes with the 32 bit ABI */
typedef int 		__CLPK_integer;
typedef int 		__CLPK_logical;
typedef float 		__CLPK_real;
typedef double 		__CLPK_doublereal;
typedef 		__CLPK_logical 	(*__CLPK_L_fp)();
typedef int 		__CLPK_ftnlen;
#else
typedef long int 	__CLPK_integer;
typedef long int 	__CLPK_logical;
typedef float 		__CLPK_real;
typedef double 		__CLPK_doublereal;
typedef __CLPK_logical 	(*__CLPK_L_fp)();
typedef long int 	__CLPK_ftnlen;
#endif

typedef struct { __CLPK_real r, i; } __CLPK_complex;
typedef struct { __CLPK_doublereal r, i; } __CLPK_doublecomplex;


#ifdef __cplusplus
}
#endif
#endif /* __CLAPACK_H */

