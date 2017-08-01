function codegen_admm(filename)

fw = fopen(filename, 'w');

% system header files
fprintf(fw, '#include <stdio.h>\n');
fprintf(fw, '#include <stdlib.h>\n');
fprintf(fw, '#include <string.h>\n');
fprintf(fw, '#include <math.h>\n');
fprintf(fw, '#include <time.h>\n\n');

% custom header files
fprintf(fw, '// custom header files\n');
fprintf(fw, '#include "projections.h"\n');
fprintf(fw, '#include "basic_operations.h"\n');
fprintf(fw, '#include "common_functions.h"\n\n');

fprintf(fw, '#include "ADMM_matrices.h"\n\n');

% main function - codegen : ADMM

fprintf(fw, '/******************************************\n');
fprintf(fw, '*          RUN ADMM ALGORITHM             *\n');
fprintf(fw, '*******************************************/\n\n');

fprintf(fw, 'int main(void)\n');
fprintf(fw, '{\n');

fprintf(fw, '\t int i;\n\n');

fprintf(fw, '\t clock_t begin, end;\n');
fprintf(fw, '\t double elapsed_time;\n\n');

fprintf(fw, '\t // temporary array to store the solution of the KKT system\n');
fprintf(fw, '\t float t[CONST_VEC_LEN]; \n\n');

fprintf(fw, '\t // primal residual\n');
fprintf(fw, '\t float rPrimal; \n\n');

fprintf(fw, '\t // dual residual\n');
fprintf(fw, '\t float rDual; \n\n');

fprintf(fw, '\t // start measure runtime \n');
fprintf(fw, '\t begin = clock();\n\n');

fprintf(fw, '\t for(i = 0;i < MAXITER; i++)\n');
fprintf(fw, '\t {\n');

fprintf(fw, '\t\t // step 1: solve linear system \n');
fprintf(fw, '\t\t solve_KKT(t); \n');
fprintf(fw, '\t\t select_values(x, t, 0, NUM_OPT_VAR);\n\n');

fprintf(fw, '\t\t // copy y into prev_y \n');
fprintf(fw, '\t\t copy_values(prev_y, y, Y_LEN); \n\n');

fprintf(fw, '\t\t // step 2: projection\n');
fprintf(fw, '\t\t compute_vector_to_project();\n');
fprintf(fw, '\t\t compute_projection();\n\n');

fprintf(fw, '\t\t // step 3: dual_update\n');
fprintf(fw, '\t\t subtract_vectors(lam, y, q, LAM_LEN);\n');
fprintf(fw, '\t\t scalar_vector_product(lam, STEPSIZE, LAM_LEN);\n\n'); 

fprintf(fw, '\t\t // step 4: convergence check \n');
fprintf(fw, '\t\t rPrimal = compute_primal_residual();\n');
fprintf(fw, '\t\t rDual   = compute_dual_residual();\n\n');

fprintf(fw, '\t\t // if solution is within the prescribed tolerance: break \n');
fprintf(fw, '\t\t if( rDual < DUALTOL && rPrimal < PRIMALTOL)\n');
fprintf(fw, '\t\t\t break;\n');
fprintf(fw, '\t } \n\n');

fprintf(fw, '\t // calculate runtime \n');
fprintf(fw, '\t end = clock(); \n');
fprintf(fw, '\t elapsed_time = (double)(end - begin) / CLOCKS_PER_SEC;\n\n');

str = sprintf('%s', 'printf("Runtime: %f ms\n", 1000*elapsed_time);');
fprintf(fw, '\t %s \n\n', str);

fprintf(fw, '\t print_the_vector(x, NUM_OPT_VAR); \n\n');
fprintf(fw, '\t return 0;\n');
fprintf(fw, '}\n\n');

fclose(fw);