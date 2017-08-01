function code_generation_FGM(type)

%% Generate the corresponding C code, particular alogrithm: ADMM
filename = 'FGM.c';
fw = fopen(filename, 'w');

%% definition of system and user defined header files + main function prototypes

fprintf(fw, '#include <stdio.h>\n');
fprintf(fw, '#include <time.h>\n');
fprintf(fw, '#include <math.h> // sqrt function \n\n');

fprintf(fw, '#include "FGM_params.h"\n\n');

fprintf(fw, '//function prototypes - main functions\n');
fprintf(fw, 'void compute_projection(void);\n');
fprintf(fw, 'void compute_roots(float *roots, float alpha_i);\n');
fprintf(fw, 'void compare_roots(float *roots, float *alpha_i_1);\n');
fprintf(fw, 'void compute_beta(float *beta_i, float alpha_i, float alpha_i_1);\n');
fprintf(fw, 'void compute_u_hat(float beta_i);\n\n');


%% auxiliary functions

fprintf(fw, '//auxiliary functions\n');
fprintf(fw, 'void add_vectors(float *resulting_vector, float *temp_vector1, float *temp_vector2, int size);\n'); 
fprintf(fw, 'void subtract_vectors(float *resulting_vector, float *temp_vector1, float *temp_vector2, int size);\n');
fprintf(fw, 'void find_minimum_element(int size);\n');
fprintf(fw, 'void find_maximal_element(int size);\n');
fprintf(fw, 'void copy_values(float *array1, float *array2, int size);\n');
fprintf(fw, 'void print_the_vector(float *A, int size);\n\n');

fprintf(fw, '//find a better a solution - i.e. redefine the pointer to a function maybe to an array of pointers\n');
fprintf(fw, 'void matrix_vector_product1(float (*A)[HG_ROWS], float *b, float *C, int size);\n');
fprintf(fw, 'void matrix_vector_product2(float (*A)[M1G_COLS], float *b, float *C, int size1, int size2);\n');
fprintf(fw, 'void scalar_vector_product(float *resulting_vector, float *vector, float scalar, int size);\n\n');


%% ADMM header 

fprintf(fw, '/******************************************\n');
fprintf(fw, '*           RUN Fast Gradient             *\n');
fprintf(fw, '*******************************************/\n\n');

%% definition of main function - variable declaration and invere computation

fprintf(fw, 'int main(void)\n');
fprintf(fw, '{\n');
fprintf(fw, '\t float x[2]; // store the roots of the quadratic equation;\n');
fprintf(fw, '\t float alpha_i_1, beta_i;\n');
fprintf(fw, '\t float alpha_i = ALPHA_I;\n');
fprintf(fw, '\t int i;\n\n');

%% main function - main loop

fprintf(fw, '\t for(i = 0;i < STEP; i++)\n');
fprintf(fw, '\t\t {\n');


fprintf(fw, '\t\t //main procedures to obtain the solution \n');
fprintf(fw, '\t\t compute_projection(); \n');
fprintf(fw, '\t\t compute_roots(x, alpha_i);\n');
fprintf(fw, '\t\t compare_roots(x, &alpha_i_1);\n');
fprintf(fw, '\t\t compute_beta(&beta_i, alpha_i, alpha_i_1);\n');
fprintf(fw, '\t\t compute_u_hat(beta_i);\n\n');

fprintf(fw, '\t\t //reassign values of alpha - convergence check\n');
fprintf(fw, '\t\t alpha_i = alpha_i_1;\n\n');

fprintf(fw, '\t\t //copy actual solution\n');
fprintf(fw, '\t\t copy_values(u_fg_prev, u_fg, HG_ROWS); \n');
fprintf(fw, '\t\t copy_values(u_fg_hat_prev, u_fg_hat, HG_ROWS);\n');

fprintf(fw, '\t\t }\n\n');

fprintf(fw, '\t //print the solution on the screen\n');
fprintf(fw, '\t print_the_vector(u_fg, HG_ROWS);\n\n');

fprintf(fw, '\t return 0;\n');
fprintf(fw, '}\n\n');

%% compute projections - functions header

fprintf(fw, '/******************************************************\n');
fprintf(fw, '*                 compute projections                 *\n');
fprintf(fw, '*******************************************************/\n\n');

%% compute projections - function definition

fprintf(fw, 'void compute_projection(void)\n');
fprintf(fw, '\t { \n');
fprintf(fw, '\t float temp_vector1[HG_ROWS]; \n');
fprintf(fw, '\t float temp_vector2[HG_ROWS]; \n');
fprintf(fw, '\t float temp_vector3[HG_ROWS]; \n');
fprintf(fw, '\t float temp_vector4[HG_ROWS]; \n\n');

fprintf(fw, '\t // compute u_fg - see matlab code: later on reformuate it \n');
fprintf(fw, '\t matrix_vector_product1(H_g, u_fg_hat_prev, temp_vector1, HG_ROWS); \n');
fprintf(fw, '\t matrix_vector_product2(M1_g, x0, temp_vector2, M1G_ROWS, M1G_COLS); \n');
fprintf(fw, '\t add_vectors(temp_vector3, temp_vector1, temp_vector2, HG_ROWS); \n');
fprintf(fw, '\t scalar_vector_product(temp_vector4, temp_vector3, 1/L, HG_ROWS); \n');
fprintf(fw, '\t subtract_vectors(u_fg, u_fg_hat_prev, temp_vector4, HG_ROWS); \n\n');

fprintf(fw, '\t // projection itself \n');
fprintf(fw, '\t find_minimum_element(HG_ROWS); \n');
fprintf(fw, '\t find_maximal_element(HG_ROWS); \n');

fprintf(fw, '} \n\n');

%% compute roots - function header

fprintf(fw, '/*******************************************************\n');
fprintf(fw, '*                     compute roots                    *\n');
fprintf(fw, '********************************************************/\n\n');

%% compute roots - function definition

fprintf(fw, 'void compute_roots(float *roots, float alpha_i) \n');
fprintf(fw, '{\n');

fprintf(fw, '\t // params in the quadratic equation \n');
fprintf(fw, '\t float a = 1; \n');
fprintf(fw, '\t float b = alpha_i * alpha_i - MU/L; \n');
fprintf(fw, '\t float c = -alpha_i * alpha_i; \n');
fprintf(fw, '\t float D = b * b - 4 * a * c; \n\n');

fprintf(fw, '\t // first root \n');
fprintf(fw, '\t roots[0] = (-1*b - sqrt(D))/2*a; \n\n');

fprintf(fw, '\t // second root \n');
fprintf(fw, '\t roots[1] = (-1*b + sqrt(D))/2*a; \n\n');

fprintf(fw, '}\n\n');

%% compare roots - function header

fprintf(fw, '/**************************************************************\n');
fprintf(fw, '*                       compare roots                         *\n');
fprintf(fw, '***************************************************************/\n\n');

%% compare roots - function definition

fprintf(fw, 'void compare_roots(float *roots, float *alpha_i_1)\n');
fprintf(fw, '{\n');

fprintf(fw, '\t if( roots[0] < 1 && roots[0] > 0) \n');
fprintf(fw, '\t\t *alpha_i_1 = roots[0]; \n');
fprintf(fw, '\t else \n');
fprintf(fw, '\t\t *alpha_i_1 = roots[1]; \n');

fprintf(fw, '}\n');

%% compute beta - function header

fprintf(fw, '/**************************************************************\n');
fprintf(fw, '*                         compute beta                        *\n');
fprintf(fw, '***************************************************************/\n\n');

%% compute beta - function definition

fprintf(fw, 'void compute_beta(float *beta_i, float alpha_i, float alpha_i_1) \n');
fprintf(fw, '{ \n');
fprintf(fw, '\t *beta_i = alpha_i * (1 - alpha_i) / (alpha_i * alpha_i + alpha_i_1); \n');
fprintf(fw, '} \n');

%% compute u_hat - function header

fprintf(fw, '/******************************************************\n');
fprintf(fw, '*                     compute u_hat                   *\n');
fprintf(fw, '*******************************************************/\n\n');

%% compute u_hat - function definition

fprintf(fw, 'void compute_u_hat(float beta_i) \n');
fprintf(fw, '{\n');

fprintf(fw, '\t float temp_vector[HG_ROWS]; \n');
fprintf(fw, '\t float temp_vector2[HG_ROWS]; \n');
fprintf(fw, '\t float beta_copy = beta_i; \n\n');

fprintf(fw, '\t subtract_vectors(temp_vector, u_fg, u_fg_prev, HG_ROWS); \n');
fprintf(fw, '\t scalar_vector_product(temp_vector2, temp_vector, beta_copy, HG_ROWS); \n');
fprintf(fw, '\t add_vectors(u_fg_hat, u_fg, temp_vector2, HG_ROWS); \n');

fprintf(fw, '}\n\n');

%% vector addition - function header

fprintf(fw, '/******************************************************\n');
fprintf(fw, '*          add two vectors of arbitrary size          *\n');
fprintf(fw, '*******************************************************/\n\n');

%% vector addition - function definition

fprintf(fw, 'void add_vectors(float *resulting_vector, float *temp_vector1, float *temp_vector2, int size)\n');
fprintf(fw, '{\n');

fprintf(fw, '\t int i;\n');
fprintf(fw, '\t for(i = 0;i < size; i++) \n');
fprintf(fw, '\t {\n');
fprintf(fw, '\t\t resulting_vector[i] = temp_vector1[i] + temp_vector2[i];\n');
fprintf(fw, '\t }\n');

fprintf(fw, '}\n\n');

%% vector subtraction - function header

fprintf(fw, '/******************************************************\n');
fprintf(fw, '*     subtract two vectors of arbitrary size          *\n');
fprintf(fw, '*******************************************************/\n\n');

%% vector subtraction - function definition

fprintf(fw, 'void subtract_vectors(float *resulting_vector, float *temp_vector1, float *temp_vector2, int size)\n');
fprintf(fw, '{\n');

fprintf(fw, '\t int i;\n');
fprintf(fw, '\t for(i = 0;i < size; i++) \n');
fprintf(fw, '\t {\n');
fprintf(fw, '\t\t resulting_vector[i] = temp_vector1[i] - temp_vector2[i];\n');
fprintf(fw, '\t }\n');

fprintf(fw, '}\n\n');

%% find minimal element - function header

fprintf(fw, '/******************************************************\n');
fprintf(fw, '* find element with minimal value in the given array  *\n');
fprintf(fw, '*******************************************************/\n\n');

%% function declaration - find minimal element

fprintf(fw, 'void find_minimum_element(int size)\n');
fprintf(fw, '{\n');

fprintf(fw, '\t int i; \n');
fprintf(fw, '\t for (i = 0; i < size; i++)\n');
fprintf(fw, '\t {\n');
fprintf(fw, '\t\t if (umax[i] <= u_fg[i])\n');
fprintf(fw, '\t\t {\n');
fprintf(fw, '\t\t\t u_fg[i] = umax[i];\n');
fprintf(fw, '\t\t }\n');
fprintf(fw, '\t }\n');

fprintf(fw, '}\n\n');

%% find maximal element - function header

fprintf(fw, '/******************************************************\n');
fprintf(fw, '* find element with maximum value in the given array  *\n');
fprintf(fw, '*******************************************************/\n\n');

%% find maximal element - function definition

fprintf(fw, 'void find_maximal_element(int size)\n');
fprintf(fw, '{\n');

fprintf(fw, '\t int i; \n');
fprintf(fw, '\t for (i = 0; i < size; i++)\n');
fprintf(fw, '\t {\n');
fprintf(fw, '\t\t if (umin[i] >= u_fg[i])\n');
fprintf(fw, '\t\t {\n');
fprintf(fw, '\t\t\t u_fg[i] = umin[i];\n');
fprintf(fw, '\t\t }\n');
fprintf(fw, '\t }\n');

fprintf(fw, '}\n\n');

%% multiply vector with a scalar - function header

fprintf(fw, '/******************************************************\n');
fprintf(fw, '*              multiply a vector with a scalar        *\n');
fprintf(fw, '*******************************************************/\n\n');

%% multiply a vector with a scalar

fprintf(fw, 'void scalar_vector_product(float *resulting_vector, float *vector, float scalar, int size) \n');
fprintf(fw, '{\n');

fprintf(fw, '\t int i; \n');
fprintf(fw, '\t for(i = 0;i < size;i++) \n');
fprintf(fw, '\t\t { \n');
fprintf(fw, '\t\t resulting_vector[i] = vector[i] * scalar; \n');
fprintf(fw, '\t\t }\n');
fprintf(fw, '}\n\n');

%% compute matrix vector product 1 - function header

fprintf(fw, '/******************************************************\n');
fprintf(fw, '*           Compute the matrix vector product         *\n');
fprintf(fw, '*******************************************************/\n\n');

%% compute matrix vector product - function definition

fprintf(fw, 'void matrix_vector_product1(float (*A)[HG_ROWS], float *b, float *C, int size) \n');
fprintf(fw, '{ \n');

fprintf(fw, '\t int i, j; \n');
fprintf(fw, '\t float sum; \n\n');

fprintf(fw, '\t for(i = 0; i < size; i++) \n');
fprintf(fw, '\t { \n');
fprintf(fw, '\t\t sum = 0.0; \n');
fprintf(fw, '\t\t for(j = 0; j < size; j++) \n');
fprintf(fw, '\t\t {\n');
fprintf(fw, '\t\t\t sum = sum + A[i][j] * b[j]; \n');
fprintf(fw, '\t\t }\n');
fprintf(fw, '\t\t C[i] = sum;\n');

fprintf(fw, '\t } \n');
fprintf(fw, '}\n\n');

%% compute matrix vector product 2 - function header

fprintf(fw, '/******************************************************\n');
fprintf(fw, '*           Compute the matrix vector product         *\n');
fprintf(fw, '*******************************************************/\n\n');

%% compute matrix vector product - function definition

fprintf(fw, 'void matrix_vector_product2(float (*A)[M1G_COLS], float *b, float *C, int size1, int size2) \n');
fprintf(fw, '{ \n');

fprintf(fw, '\t int i, j; \n');
fprintf(fw, '\t float sum; \n\n');

fprintf(fw, '\t for(i = 0; i < size1; i++) \n');
fprintf(fw, '\t { \n');
fprintf(fw, '\t\t sum = 0.0; \n');
fprintf(fw, '\t\t for(j = 0; j < size2; j++) \n');
fprintf(fw, '\t\t {\n');
fprintf(fw, '\t\t\t sum = sum + A[i][j] * b[j]; \n');
fprintf(fw, '\t\t }\n');
fprintf(fw, '\t\t C[i] = sum;\n');

fprintf(fw, '\t } \n');
fprintf(fw, '}\n\n');

%% copy values - function header

fprintf(fw, '/******************************************************\n');
fprintf(fw, '*   save the previous values of the opt. variables    *\n');
fprintf(fw, '*******************************************************/\n\n');

%% function declaration - copy values

fprintf(fw, 'void copy_values(float *array1, float *array2, int size)\n');
fprintf(fw, '{\n');

fprintf(fw, '\t int i; \n');
fprintf(fw, '\t for(i = 0;i < size;i++) \n');
fprintf(fw, '\t\t { \n');
fprintf(fw, '\t\t array1[i] = array2[i]; \n');
fprintf(fw, '\t\t }\n');
fprintf(fw, '}\n\n');

%% print vector to the screen - function header

fprintf(fw, '/******************************************************\n');
fprintf(fw, '*           print the vector to the screen            *\n');
fprintf(fw, '*******************************************************/\n\n');

%% function definition - print vector to the screen

fprintf(fw, 'void print_the_vector(float *vector, int size)\n');
fprintf(fw, '{ \n');
fprintf(fw, '\t int i; \n');
fprintf(fw, '\t for(i = 0; i < size;i++)\n');
fprintf(fw, '\t { \n');
str = 'printf("%f\n", vector[i]);';
fprintf(fw, '\t\t %s \n', str);
fprintf(fw, '\t } \n');
fprintf(fw, '} \n');