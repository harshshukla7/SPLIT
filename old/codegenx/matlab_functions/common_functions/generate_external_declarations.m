function generate_external_declarations(data_ADMM)

% create a file handler and the corresponding header file
filename = 'external_declarations.h';
fw = fopen(filename, 'w');

% do not allow mulitple includes
fprintf(fw, '#ifndef _EXTERNAL_DECLARATIONS__H__\n');
fprintf(fw, '#define _EXTERNAL_DECLARATIONS__H__\n\n');

% macro ADMM: I need it for proper codegen
fprintf(fw,'#define ADMM\n');

% basic optimization data
fprintf(fw,'extern int NUM_OPT_VAR; ');
fprintf(fw,'// number of optimisation variables\n');
fprintf(fw,'extern int NUM_OF_STATES; ');
fprintf(fw,'// number of states\n');
fprintf(fw,'extern int NPROX; ');
fprintf(fw,'// number of proxes\n');
fprintf(fw,'extern float DUALTOL; ');
fprintf(fw,'// dual residual error\n');
fprintf(fw,'extern float PRIMALTOL; ');
fprintf(fw,'// primal residual error\n');
fprintf(fw,'extern int MAXITER; ');
fprintf(fw,'// maximum number of iterations\n');
fprintf(fw,'extern int TERM_CHECK_FREQ; ');
fprintf(fw,'//  frequency of termination check\n');
fprintf(fw,'extern int STEPSIZE; ');
fprintf(fw,'//  frequency of termination check\n');

% length of vectors
fprintf(fw,'extern int B_LEN; ');
fprintf(fw,'// size of the vector b\n');
fprintf(fw,'extern int F_LEN; ');
fprintf(fw,'// size of the vector f\n'); 
fprintf(fw,'extern int C_LEN; ');
fprintf(fw,'// size of the vector c\n');
fprintf(fw,'extern int L_LEN; ');
fprintf(fw,'// size of the vector l\n');
fprintf(fw,'extern int X_LEN; ');
fprintf(fw,'// size of the vector x\n');
fprintf(fw,'extern int Y_LEN; ');
fprintf(fw,'// size of the vector y\n');
fprintf(fw,'extern int LAM_LEN; ');
fprintf(fw,'// size of the vector lam\n');
fprintf(fw,'extern int Q_LEN; ');
fprintf(fw,'// size of the vector q\n');
fprintf(fw,'extern int CONST_VEC_LEN; ');
fprintf(fw,'// size of the constant vector \n');

% sizes of matrices
fprintf(fw,'extern int A_ROWS; ',                   data_ADMM.A_ROWS);
fprintf(fw,'// number of rows in matrix A\n');
fprintf(fw,'extern int A_COLS; ',                   data_ADMM.A_COLS);
fprintf(fw,'// number of columns in matrix A\n');
fprintf(fw,'extern int Q_ROWS; ',                   data_ADMM.Q_ROWS);
fprintf(fw,'// number of rows in matrix Q\n');
fprintf(fw,'extern int Q_COLS; ',                   data_ADMM.Q_COLS);
fprintf(fw,'// number of columns matrix Q\n');
fprintf(fw,'extern int L_ROWS; ',                   data_ADMM.L_ROWS);
fprintf(fw,'// number of rows in matrix L\n');
fprintf(fw,'extern int L_COLS; ',                   data_ADMM.L_COLS);
fprintf(fw,'// number of columns in matrix L\n');
fprintf(fw,'extern int SCALEDL_ROWS ; ',            data_ADMM.scaledL_ROWS);
fprintf(fw,'// number of rows in matrix rho*L\n');
fprintf(fw,'extern int SCALEDL_COLS; ',             data_ADMM.scaledL_COLS);
fprintf(fw,'// number of cols in matrix rho*L\n');
fprintf(fw,'extern int MATRIX_ROWS; ',              data_ADMM.matrix_ROWS);
fprintf(fw,'// number of rows in matrix KKT\n');
fprintf(fw,'extern int MATRIX_COLS; ',              data_ADMM.matrix_COLS);
fprintf(fw,'// number of cols in matrix KKT\n');

% sizes of matrices, containing indices of nonzero elements
fprintf(fw,'extern int  AI_ROWS; ',                  data_ADMM.AI_ROWS);
fprintf(fw,'// number of rows in matrix AI\n');
fprintf(fw,'extern int  AI_COLS; ',                  data_ADMM.AI_COLS);
fprintf(fw,'// number of columns in matrix AI\n');
fprintf(fw,'extern int  QI_ROWS; ',                  data_ADMM.QI_ROWS);
fprintf(fw,'// number of rows in matrix QI\n');
fprintf(fw,'extern int  QI_COLS; ',                  data_ADMM.QI_COLS);
fprintf(fw,'// number of columns matrix QI\n');
fprintf(fw,'extern int  LI_ROWS; ',                  data_ADMM.LI_ROWS);
fprintf(fw,'// number of rows in matrix LI\n');
fprintf(fw,'extern int  LI_COLS; ',                  data_ADMM.LI_COLS);
fprintf(fw,'// number of columns in matrix LI\n');
fprintf(fw,'extern int MATRIXI_ROWS; ',              data_ADMM.matrixI_ROWS);
fprintf(fw,'// number of rows in matrix KKT\n');
fprintf(fw,'extern int MATRIXI_COLS ; ',             data_ADMM.matrixI_COLS);
fprintf(fw,'// number of cols in matrix KKT\n');
fprintf(fw,'extern int SCALEDLI_ROWS; ',             data_ADMM.scaledLI_ROWS);
fprintf(fw,'// number of rows in matrix rho*LI\n');
fprintf(fw,'extern int SCALEDLI_COLS; ',             data_ADMM.scaledLI_COLS);
fprintf(fw,'// number of cols in matrix rho*LI\n\n');

% export 1D arrays to a header file
generate_external_variables_vectors(filename, 'x',              'extern float',  data_ADMM.x_LEN);
generate_external_variables_vectors(filename, 'y',              'extern float',  data_ADMM.y_LEN);
generate_external_variables_vectors(filename, 'prev_y',         'extern float',  data_ADMM.y_LEN);
generate_external_variables_vectors(filename, 'lam',            'extern float',  data_ADMM.lam_LEN);
generate_external_variables_vectors(filename, 'q',              'extern float',  data_ADMM.q_LEN);
generate_external_variables_vectors(filename, 'b',              'extern float',  data_ADMM.b_LEN);
generate_external_variables_vectors(filename, 'f',              'extern float',  data_ADMM.f_LEN);
generate_external_variables_vectors(filename, 'c',              'extern float',  data_ADMM.c_LEN);
generate_external_variables_vectors(filename, 'l',              'extern float',  data_ADMM.l_LEN);
generate_external_variables_vectors(filename, 'const_vec',      'extern float',  data_ADMM.const_vec_LEN);
generate_external_variables_vectors(filename, 'new_prox_beg',   'extern float',  data_ADMM.nProx);
generate_external_variables_vectors(filename, 'new_prox_end',   'extern float',  data_ADMM.nProx);
generate_external_variables_vectors(filename, 'len_of_vectors', 'extern float',  data_ADMM.nProx);

% export matrices to a header file
generate_external_variables_matrices(filename, 'A',        'extern float',  data_ADMM.A_ROWS,       data_ADMM.A_COLS);
generate_external_variables_matrices(filename, 'Q',        'extern float',  data_ADMM.Q_ROWS,       data_ADMM.Q_COLS);
generate_external_variables_matrices(filename, 'L',        'extern float',  data_ADMM.L_ROWS,       data_ADMM.L_COLS);
generate_external_variables_matrices(filename, 'matrix',   'extern float',  data_ADMM.matrix_ROWS,  data_ADMM.matrix_COLS);
generate_external_variables_matrices(filename, 'scaledL',  'extern float',  data_ADMM.scaledL_ROWS, data_ADMM.scaledL_COLS);

% export matrice to a header file - indices of nonzero elements
generate_external_variables_matrices(filename, 'AI',        'extern int',  data_ADMM.AI_ROWS,        data_ADMM.AI_COLS);
generate_external_variables_matrices(filename, 'QI',        'extern int',  data_ADMM.QI_ROWS,        data_ADMM.QI_COLS);
generate_external_variables_matrices(filename, 'LI',        'extern int',  data_ADMM.LI_ROWS,        data_ADMM.LI_COLS);
generate_external_variables_matrices(filename, 'matrixI',   'extern int',  data_ADMM.matrixI_ROWS,   data_ADMM.matrixI_COLS);
generate_external_variables_matrices(filename, 'scaledLI',  'extern int',  data_ADMM.scaledLI_ROWS,  data_ADMM.scaledLI_COLS);

% generate array of proxes and norms
generate_external_array_of_prox_types(filename);
generate_external_array_of_norm_types(filename);

% generate 1D arrays containing normTypesm, weights and constants to perform projections
generate_external_variables_vectors(filename,   'weights',   'extern float',  data_ADMM.nProx);
generate_external_variables_vectors(filename, 'constants',   'extern float',  data_ADMM.nProx);

% preprocessor directive ends
fw = fopen(filename, 'a');
fprintf(fw, '\n');
fprintf(fw, '#endif');

% close the file handler, otherwise std::exception is thrown
fclose(fw);