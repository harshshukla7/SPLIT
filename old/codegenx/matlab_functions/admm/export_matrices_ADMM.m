function export_matrices_ADMM(data_ADMM)

% create a file handler and the corresponding header file
filename = 'ADMM_matrices.h';
fw = fopen(filename, 'w');

% do not allow mulitple includes
fprintf(fw, '#ifndef _ADMM_MATRICES__H__\n');
fprintf(fw, '#define _ADMM_MATRICES__H__\n\n');

% basic optimization data
fprintf(fw,'int NUM_OPT_VAR = %d; ',                data_ADMM.num_of_opt_var);
fprintf(fw,'// number of optimisation variables\n');
fprintf(fw,'int NUM_OF_STATES = %d; ',              data_ADMM.num_of_states);
fprintf(fw,'// number of states\n');
fprintf(fw,'int NPROX  = %d; ',                     data_ADMM.nProx);
fprintf(fw,'// number of proxes\n');
fprintf(fw,'double DUALTOL =  %.15f; ',             data_ADMM.settings.dualTol);
fprintf(fw,'// dual residual error\n');
fprintf(fw,'double PRIMALTOL = %.15f; ',            data_ADMM.settings.primalTol);
fprintf(fw,'// primal residual error\n');
fprintf(fw,'double ALPHA = %.15f; ',                data_ADMM.settings.relaxation);
fprintf(fw,'// parameter for over relaxation\n');
fprintf(fw,'int MAXITER  =  %d; ',                  data_ADMM.settings.maxItr);
fprintf(fw,'// maximum number of iterations\n');
fprintf(fw,'int TERM_CHECK_FREQ =  %d; ',           data_ADMM.settings.terminationCheckFreq);
fprintf(fw,'// frequency of termination check\n');
fprintf(fw,'double STEPSIZE = %.15f; ',             data_ADMM.settings.rho);
fprintf(fw,'// frequency of termination check\n');

% length of vectors
fprintf(fw,'int B_LEN  =  %d; ',                    data_ADMM.b_LEN);
fprintf(fw,'// size of the vector b\n');
fprintf(fw,'int F_LEN =   %d; ',                    data_ADMM.f_LEN);
fprintf(fw,'// size of the vector f\n'); 
fprintf(fw,'int C_LEN =  %d; ',                     data_ADMM.c_LEN);
fprintf(fw,'// size of the vector c\n');
fprintf(fw,'int L_LEN = %d; ',                      data_ADMM.l_LEN);
fprintf(fw,'// size of the vector l\n');
fprintf(fw,'int X_LEN = %d; ',                      data_ADMM.x_LEN);
fprintf(fw,'// size of the vector x\n');
fprintf(fw,'int Y_LEN = %d; ',                      data_ADMM.y_LEN);
fprintf(fw,'// size of the vector y\n');
fprintf(fw,'int LAM_LEN = %d; ',                    data_ADMM.lam_LEN);
fprintf(fw,'// size of the vector lam\n');
fprintf(fw,'int Q_LEN =  %d; ',                     data_ADMM.q_LEN);
fprintf(fw,'// size of the vector q\n');
fprintf(fw,'int T_LEN =  %d; ',                     data_ADMM.t_LEN);
fprintf(fw,'// size of the vector t\n');
fprintf(fw,'int CONST_VEC_LEN =  %d; ',             data_ADMM.const_vec_LEN);
fprintf(fw,'// size of the constant vector \n');
fprintf(fw,'int LX_LEN =  %d; ',                    data_ADMM.Lx_LEN);
fprintf(fw,'// storage for over relaxation \n');

% sizes of matrices
fprintf(fw,'int A_ROWS = %d; ',                     data_ADMM.A_ROWS);
fprintf(fw,'// number of rows in matrix A\n');
fprintf(fw,'int A_COLS = %d; ',                     data_ADMM.A_COLS);
fprintf(fw,'// number of columns in matrix A\n');
fprintf(fw,'int Q_ROWS = %d; ',                     data_ADMM.Q_ROWS);
fprintf(fw,'// number of rows in matrix Q\n');
fprintf(fw,'int Q_COLS = %d; ',                     data_ADMM.Q_COLS);
fprintf(fw,'// number of columns matrix Q\n');
fprintf(fw,'int L_ROWS = %d; ',                     data_ADMM.L_ROWS);
fprintf(fw,'// number of rows in matrix L\n');
fprintf(fw,'int L_COLS = %d; ',                     data_ADMM.L_COLS);
fprintf(fw,'// number of columns in matrix L\n');
fprintf(fw,'int SCALEDL_ROWS = %d; ',               data_ADMM.scaledL_ROWS);
fprintf(fw,'// number of rows in matrix rho*L\n');
fprintf(fw,'int SCALEDL_COLS = %d; ',               data_ADMM.scaledL_COLS);
fprintf(fw,'// number of cols in matrix rho*L\n');
fprintf(fw,'int MATRIX_ROWS =  %d; ',               data_ADMM.matrix_ROWS);
fprintf(fw,'// number of rows in matrix KKT\n');
fprintf(fw,'int MATRIX_COLS = %d; ',                data_ADMM.matrix_COLS);
fprintf(fw,'// number of cols in matrix KKT\n');
fprintf(fw,'int ALPHAL_ROWS = %d; ',                data_ADMM.alphaL_ROWS);
fprintf(fw,'// number of rows in matrix alpha*L\n');
fprintf(fw,'int ALPHAL_COLS = %d; ',                data_ADMM.alphaL_COLS);
fprintf(fw,'// number of cols in matrix alpha*L\n');

% sizes of matrices, containing indices of nonzero elements
fprintf(fw,'int  AI_ROWS =  %d; ',                  data_ADMM.AI_ROWS);
fprintf(fw,'// number of rows in matrix AI\n');
fprintf(fw,'int  AI_COLS =  %d; ',                  data_ADMM.AI_COLS);
fprintf(fw,'// number of columns in matrix AI\n');
fprintf(fw,'int  QI_ROWS =  %d; ',                  data_ADMM.QI_ROWS);
fprintf(fw,'// number of rows in matrix QI\n');
fprintf(fw,'int  QI_COLS =  %d; ',                  data_ADMM.QI_COLS);
fprintf(fw,'// number of columns matrix QI\n');
fprintf(fw,'int  LI_ROWS =  %d; ',                  data_ADMM.LI_ROWS);
fprintf(fw,'// number of rows in matrix LI\n');
fprintf(fw,'int  LI_COLS =  %d; ',                  data_ADMM.LI_COLS);
fprintf(fw,'// number of columns in matrix LI\n');
fprintf(fw,'int MATRIXI_ROWS =  %d; ',              data_ADMM.matrixI_ROWS);
fprintf(fw,'// number of rows in matrix KKT\n');
fprintf(fw,'int MATRIXI_COLS =  %d; ',              data_ADMM.matrixI_COLS);
fprintf(fw,'// number of cols in matrix KKT\n');
fprintf(fw,'int SCALEDLI_ROWS =  %d; ',             data_ADMM.scaledLI_ROWS);
fprintf(fw,'// number of rows in matrix rho*LI\n');
fprintf(fw,'int SCALEDLI_COLS =  %d; ',             data_ADMM.scaledLI_COLS);
fprintf(fw,'// number of cols in matrix rho*LI\n');
fprintf(fw,'int ALPHALI_ROWS =  %d; ',              data_ADMM.alphaLI_ROWS);
fprintf(fw,'// number of rows in matrix alpha*LI\n');
fprintf(fw,'int ALPHALI_COLS =  %d; ',              data_ADMM.alphaLI_COLS);
fprintf(fw,'// number of cols in matrix alpha*LI\n\n');

% export 1D arrays to a header file
generate_vector(filename, data_ADMM.x,              'x',              'double',  data_ADMM.x_LEN);
generate_vector(filename, data_ADMM.y,              'y',              'double',  data_ADMM.y_LEN);
generate_vector(filename, data_ADMM.prev_y,         'prev_y',         'double',  data_ADMM.y_LEN);
generate_vector(filename, data_ADMM.lam,            'lam',            'double',  data_ADMM.lam_LEN);
generate_vector(filename, data_ADMM.q,              'q',              'double',  data_ADMM.q_LEN);
generate_vector(filename, data_ADMM.b,              'b',              'double',  data_ADMM.b_LEN);
generate_vector(filename, data_ADMM.f,              'f',              'double',  data_ADMM.f_LEN);
generate_vector(filename, data_ADMM.c,              'c',              'double',  data_ADMM.c_LEN);
generate_vector(filename, data_ADMM.l,              'l',              'double',  data_ADMM.l_LEN);
generate_vector(filename, data_ADMM.const_vec,      'const_vec',      'double',  data_ADMM.const_vec_LEN);
generate_vector(filename, data_ADMM.new_prox_beg,   'new_prox_beg',   'int',     data_ADMM.nProx);
generate_vector(filename, data_ADMM.new_prox_end,   'new_prox_end',   'int',     data_ADMM.nProx);
generate_vector(filename, data_ADMM.t           ,   't',              'double',  data_ADMM.t_LEN);
generate_vector(filename, data_ADMM.lx          ,   'lx',             'double',  data_ADMM.Lx_LEN);
generate_vector(filename, data_ADMM.len_of_vectors, 'len_of_vectors', 'int',     data_ADMM.nProx);

% export matrices to a header file
generate_matrix(filename, data_ADMM.A,        'A',        'double',  data_ADMM.A_ROWS,       data_ADMM.A_COLS);
generate_matrix(filename, data_ADMM.Q,        'Q',        'double',  data_ADMM.Q_ROWS,       data_ADMM.Q_COLS);
generate_matrix(filename, data_ADMM.L,        'L',        'double',  data_ADMM.L_ROWS,       data_ADMM.L_COLS);
generate_matrix(filename, data_ADMM.matrix,   'matrix',   'double',  data_ADMM.matrix_ROWS,  data_ADMM.matrix_COLS);
generate_matrix(filename, data_ADMM.scaledL,  'scaledL',  'double',  data_ADMM.scaledL_ROWS, data_ADMM.scaledL_COLS);
generate_matrix(filename, data_ADMM.alphaL,   'alphaL',   'double',  data_ADMM.alphaL_ROWS,  data_ADMM.alphaL_COLS);

% export matrice to a header file - indices of nonzero elements
generate_matrix(filename, data_ADMM.AI,        'AI',        'int',  data_ADMM.AI_ROWS,        data_ADMM.AI_COLS);
generate_matrix(filename, data_ADMM.QI,        'QI',        'int',  data_ADMM.QI_ROWS,        data_ADMM.QI_COLS);
generate_matrix(filename, data_ADMM.LI,        'LI',        'int',  data_ADMM.LI_ROWS,        data_ADMM.LI_COLS);
generate_matrix(filename, data_ADMM.matrixI,   'matrixI',   'int',  data_ADMM.matrixI_ROWS,   data_ADMM.matrixI_COLS);
generate_matrix(filename, data_ADMM.scaledLI,  'scaledLI',  'int',  data_ADMM.scaledLI_ROWS,  data_ADMM.scaledLI_COLS);
generate_matrix(filename, data_ADMM.alphaLI,   'alphaLI',   'int',  data_ADMM.alphaLI_ROWS,   data_ADMM.alphaLI_COLS);

% generate array of proxes and norms
generate_array_of_prox_types(filename, data_ADMM.prox_type, data_ADMM.nProx);
generate_array_of_norm_types(filename, data_ADMM.normType,  data_ADMM.nProx);

% generate 1D arrays containing normTypesm, weights and constants to perform projections
generate_vector(filename, data_ADMM.weight,    'weights',   'double',  data_ADMM.nProx);
generate_vector(filename, data_ADMM.constants, 'constants', 'double',  data_ADMM.nProx);

% preprocessor directive ends
fw = fopen(filename, 'a');
fprintf(fw, '\n');
fprintf(fw, '#endif');

% close the file handler, otherwise std::exception is thrown
fclose(fw);