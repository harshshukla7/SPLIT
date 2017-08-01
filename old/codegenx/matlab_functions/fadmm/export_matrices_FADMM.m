function export_matrices_FADMM(data_FADMM)

% create a file handler and the corresponding header file
filename = 'FADMM_matrices.h';
fw = fopen(filename, 'w');

% do not allow mulitple includes
fprintf(fw, '#ifndef _FADMM_MATRICES__H__\n');
fprintf(fw, '#define _FADMM_MATRICES__H__\n\n');

% basic optimization data
fprintf(fw,'int NUM_OPT_VAR = %d; ',                data_FADMM.num_of_opt_var);
fprintf(fw,'// number of optimisation variables\n');
fprintf(fw,'int NUM_OF_STATES = %d; ',              data_FADMM.num_of_states);
fprintf(fw,'// number of states\n');
fprintf(fw,'int NPROX  = %d; ',                     data_FADMM.nProx);
fprintf(fw,'// number of proxes\n');
fprintf(fw,'double DUALTOL =  %.15f; ',             data_FADMM.settings.dualTol);
fprintf(fw,'// dual residual error\n');
fprintf(fw,'double PRIMALTOL = %.15f; ',            data_FADMM.settings.primalTol);
fprintf(fw,'// primal residual error\n');
fprintf(fw,'double ALPHA = %.15f; ',                data_FADMM.settings.relaxation);
fprintf(fw,'// parameter for over relaxation\n');
fprintf(fw,'int MAXITER  =  %d; ',                  data_FADMM.settings.maxItr);
fprintf(fw,'// maximum number of iterations\n');
fprintf(fw,'int TERM_CHECK_FREQ =  %d; ',           data_FADMM.settings.terminationCheckFreq);
fprintf(fw,'// frequency of termination check\n');
fprintf(fw,'double STEPSIZE = %.15f; ',             data_FADMM.settings.rho);
fprintf(fw,'// stepsize\n');

% length of vectors
fprintf(fw,'int B_LEN  =  %d; ',                    data_FADMM.b_LEN);
fprintf(fw,'// size of the vector b\n');
fprintf(fw,'int F_LEN =   %d; ',                    data_FADMM.f_LEN);
fprintf(fw,'// size of the vector f\n'); 
fprintf(fw,'int C_LEN =  %d; ',                     data_FADMM.c_LEN);
fprintf(fw,'// size of the vector c\n');
fprintf(fw,'int L_LEN = %d; ',                      data_FADMM.l_LEN);
fprintf(fw,'// size of the vector l\n');
fprintf(fw,'int X_LEN = %d; ',                      data_FADMM.x_LEN);
fprintf(fw,'// size of the vector x\n');
fprintf(fw,'int Y_LEN = %d; ',                      data_FADMM.y_LEN);
fprintf(fw,'// size of the vector y\n');
fprintf(fw,'int LAM_LEN = %d; ',                    data_FADMM.lam_LEN);
fprintf(fw,'// size of the vector lam\n');
fprintf(fw,'int Q_LEN =  %d; ',                     data_FADMM.q_LEN);
fprintf(fw,'// size of the vector q\n');
fprintf(fw,'int T_LEN =  %d; ',                     data_FADMM.t_LEN);
fprintf(fw,'// size of the vector t\n');
fprintf(fw,'int CONST_VEC_LEN =  %d; ',             data_FADMM.const_vec_LEN);
fprintf(fw,'// size of the constant vector \n');
fprintf(fw,'int LX_LEN =  %d; ',                    data_FADMM.Lx_LEN);
fprintf(fw,'// storage for over relaxation \n');

% sizes of matrices
fprintf(fw,'int A_ROWS = %d; ',                     data_FADMM.A_ROWS);
fprintf(fw,'// number of rows in matrix A\n');
fprintf(fw,'int A_COLS = %d; ',                     data_FADMM.A_COLS);
fprintf(fw,'// number of columns in matrix A\n');
fprintf(fw,'int Q_ROWS = %d; ',                     data_FADMM.Q_ROWS);
fprintf(fw,'// number of rows in matrix Q\n');
fprintf(fw,'int Q_COLS = %d; ',                     data_FADMM.Q_COLS);
fprintf(fw,'// number of columns matrix Q\n');
fprintf(fw,'int L_ROWS = %d; ',                     data_FADMM.L_ROWS);
fprintf(fw,'// number of rows in matrix L\n');
fprintf(fw,'int L_COLS = %d; ',                     data_FADMM.L_COLS);
fprintf(fw,'// number of columns in matrix L\n');
fprintf(fw,'int SCALEDL_ROWS = %d; ',               data_FADMM.scaledL_ROWS);
fprintf(fw,'// number of rows in matrix rho*L\n');
fprintf(fw,'int SCALEDL_COLS = %d; ',               data_FADMM.scaledL_COLS);
fprintf(fw,'// number of cols in matrix rho*L\n');
fprintf(fw,'int MATRIX_ROWS =  %d; ',               data_FADMM.matrix_ROWS);
fprintf(fw,'// number of rows in matrix KKT\n');
fprintf(fw,'int MATRIX_COLS = %d; ',                data_FADMM.matrix_COLS);
fprintf(fw,'// number of cols in matrix KKT\n');
fprintf(fw,'int ALPHAL_ROWS = %d; ',                data_FADMM.alphaL_ROWS);
fprintf(fw,'// number of rows in matrix alpha*L\n');
fprintf(fw,'int ALPHAL_COLS = %d; ',                data_FADMM.alphaL_COLS);
fprintf(fw,'// number of cols in matrix alpha*L\n');


% sizes of matrices, containing indices of nonzero elements
fprintf(fw,'int  AI_ROWS =  %d; ',                  data_FADMM.AI_ROWS);
fprintf(fw,'// number of rows in matrix AI\n');
fprintf(fw,'int  AI_COLS =  %d; ',                  data_FADMM.AI_COLS);
fprintf(fw,'// number of columns in matrix AI\n');
fprintf(fw,'int  QI_ROWS =  %d; ',                  data_FADMM.QI_ROWS);
fprintf(fw,'// number of rows in matrix QI\n');
fprintf(fw,'int  QI_COLS =  %d; ',                  data_FADMM.QI_COLS);
fprintf(fw,'// number of columns matrix QI\n');
fprintf(fw,'int  LI_ROWS =  %d; ',                  data_FADMM.LI_ROWS);
fprintf(fw,'// number of rows in matrix LI\n');
fprintf(fw,'int  LI_COLS =  %d; ',                  data_FADMM.LI_COLS);
fprintf(fw,'// number of columns in matrix LI\n');
fprintf(fw,'int MATRIXI_ROWS =  %d; ',              data_FADMM.matrixI_ROWS);
fprintf(fw,'// number of rows in matrix KKT\n');
fprintf(fw,'int MATRIXI_COLS =  %d; ',              data_FADMM.matrixI_COLS);
fprintf(fw,'// number of cols in matrix KKT\n');
fprintf(fw,'int SCALEDLI_ROWS =  %d; ',             data_FADMM.scaledLI_ROWS);
fprintf(fw,'// number of rows in matrix rho*LI\n');
fprintf(fw,'int SCALEDLI_COLS =  %d; ',             data_FADMM.scaledLI_COLS);
fprintf(fw,'// number of cols in matrix rho*LI\n');
fprintf(fw,'int ALPHALI_ROWS =  %d; ',              data_FADMM.alphaLI_ROWS);
fprintf(fw,'// number of rows in matrix alpha*LI\n');
fprintf(fw,'int ALPHALI_COLS =  %d; ',              data_FADMM.alphaLI_COLS);
fprintf(fw,'// number of cols in matrix alpha*LI\n\n');

% export 1D arrays to a header file
generate_vector(filename, data_FADMM.x,              'x',              'double',  data_FADMM.x_LEN);
generate_vector(filename, data_FADMM.y,              'y',              'double',  data_FADMM.y_LEN);
generate_vector(filename, data_FADMM.prev_y,         'prev_y',         'double',  data_FADMM.y_LEN);
generate_vector(filename, data_FADMM.lam,            'lam',            'double',  data_FADMM.lam_LEN);
generate_vector(filename, data_FADMM.prev_lam,       'prev_lam',       'double',  data_FADMM.lam_LEN);
generate_vector(filename, data_FADMM.lam_hat,        'lam_hat',        'double',  data_FADMM.lam_LEN);
generate_vector(filename, data_FADMM.lam_hat,        'y_hat',          'double',  data_FADMM.y_LEN);
generate_vector(filename, data_FADMM.q,              'q',              'double',  data_FADMM.q_LEN);
generate_vector(filename, data_FADMM.b,              'b',              'double',  data_FADMM.b_LEN);
generate_vector(filename, data_FADMM.f,              'f',              'double',  data_FADMM.f_LEN);
generate_vector(filename, data_FADMM.c,              'c',              'double',  data_FADMM.c_LEN);
generate_vector(filename, data_FADMM.l,              'l',              'double',  data_FADMM.l_LEN);
generate_vector(filename, data_FADMM.const_vec,      'const_vec',      'double',  data_FADMM.const_vec_LEN);
generate_vector(filename, data_FADMM.new_prox_beg,   'new_prox_beg',   'int',     data_FADMM.nProx);
generate_vector(filename, data_FADMM.new_prox_end,   'new_prox_end',   'int',     data_FADMM.nProx);
generate_vector(filename, data_FADMM.t,              't',              'double',  data_FADMM.t_LEN);
generate_vector(filename, data_FADMM.lx,             'lx',             'double',  data_FADMM.Lx_LEN);
generate_vector(filename, data_FADMM.len_of_vectors, 'len_of_vectors', 'int',     data_FADMM.nProx);

% export matrices to a header file
generate_matrix(filename, data_FADMM.A,        'A',        'double',  data_FADMM.A_ROWS,       data_FADMM.A_COLS);
generate_matrix(filename, data_FADMM.Q,        'Q',        'double',  data_FADMM.Q_ROWS,       data_FADMM.Q_COLS);
generate_matrix(filename, data_FADMM.L,        'L',        'double',  data_FADMM.L_ROWS,       data_FADMM.L_COLS);
generate_matrix(filename, data_FADMM.matrix,   'matrix',   'double',  data_FADMM.matrix_ROWS,  data_FADMM.matrix_COLS);
generate_matrix(filename, data_FADMM.scaledL,  'scaledL',  'double',  data_FADMM.scaledL_ROWS, data_FADMM.scaledL_COLS);
generate_matrix(filename, data_FADMM.alphaL,   'alphaL',   'double',  data_FADMM.alphaL_ROWS,  data_FADMM.alphaL_COLS);


% export matrice to a header file - indices of nonzero elements
generate_matrix(filename, data_FADMM.AI,        'AI',        'int',  data_FADMM.AI_ROWS,        data_FADMM.AI_COLS);
generate_matrix(filename, data_FADMM.QI,        'QI',        'int',  data_FADMM.QI_ROWS,        data_FADMM.QI_COLS);
generate_matrix(filename, data_FADMM.LI,        'LI',        'int',  data_FADMM.LI_ROWS,        data_FADMM.LI_COLS);
generate_matrix(filename, data_FADMM.matrixI,   'matrixI',   'int',  data_FADMM.matrixI_ROWS,   data_FADMM.matrixI_COLS);
generate_matrix(filename, data_FADMM.scaledLI,  'scaledLI',  'int',  data_FADMM.scaledLI_ROWS,  data_FADMM.scaledLI_COLS);
generate_matrix(filename, data_FADMM.alphaLI,   'alphaLI',   'int',  data_FADMM.alphaLI_ROWS,   data_FADMM.alphaLI_COLS);

% generate array of proxes and norms
generate_array_of_prox_types(filename, data_FADMM.prox_type, data_FADMM.nProx);
generate_array_of_norm_types(filename, data_FADMM.normType,  data_FADMM.nProx);

% generate 1D arrays containing normTypesm, weights and constants to perform projections
generate_vector(filename, data_FADMM.weight,    'weights',   'double',  data_FADMM.nProx);
generate_vector(filename, data_FADMM.constants, 'constants', 'double',  data_FADMM.nProx);

% preprocessor directive ends
fw = fopen(filename, 'a');
fprintf(fw, '\n');
fprintf(fw, '#endif');

% close the file handler, otherwise std::exception is thrown
fclose(fw);