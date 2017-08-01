function export_matrices_FAMA(data_FAMA)

% create a file handler and the corresponding header file
filename = 'FAMA_matrices.h';
fw = fopen(filename, 'w');

% do not allow mulitple includes
fprintf(fw, '#ifndef _FAMA_MATRICES__H__\n');
fprintf(fw, '#define _FAMA_MATRICES__H__\n\n');

% basic optimization data
fprintf(fw,'int NUM_OPT_VAR = %d; ',                data_FAMA.num_of_opt_var);
fprintf(fw,'// number of optimisation variables\n');
fprintf(fw,'int NUM_OF_STATES = %d; ',              data_FAMA.num_of_states);
fprintf(fw,'// number of states\n');
fprintf(fw,'int NPROX  = %d; ',                     data_FAMA.nProx);
fprintf(fw,'// number of proxes\n');
fprintf(fw,'double DUALTOL =  %.15f; ',             data_FAMA.settings.dualTol);
fprintf(fw,'// dual residual error\n');
% fprintf(fw,'double PRIMALTOL = %.15f; ',            data_FAMA.settings.primalTol);
% fprintf(fw,'// primal residual error\n');
fprintf(fw,'double ALPHA = %.15f; ',                data_FAMA.settings.relaxation);
fprintf(fw,'// parameter for over relaxation\n');
fprintf(fw,'int MAXITER  =  %d; ',                  data_FAMA.settings.maxItr);
fprintf(fw,'// maximum number of iterations\n');
fprintf(fw,'int TERM_CHECK_FREQ =  %d; ',           data_FAMA.settings.terminationCheckFreq);
fprintf(fw,'// frequency of termination check\n');
fprintf(fw,'double STEPSIZE = %.15f; ',             data_FAMA.settings.rho);
fprintf(fw,'// stepsize\n');

% length of vectors
fprintf(fw,'int B_LEN  =  %d; ',                    data_FAMA.b_LEN);
fprintf(fw,'// size of the vector b\n');
fprintf(fw,'int F_LEN =   %d; ',                    data_FAMA.f_LEN);
fprintf(fw,'// size of the vector f\n'); 
fprintf(fw,'int C_LEN =  %d; ',                     data_FAMA.c_LEN);
fprintf(fw,'// size of the vector c\n');
fprintf(fw,'int L_LEN = %d; ',                      data_FAMA.l_LEN);
fprintf(fw,'// size of the vector l\n');
fprintf(fw,'int X_LEN = %d; ',                      data_FAMA.x_LEN);
fprintf(fw,'// size of the vector x\n');
fprintf(fw,'int Y_LEN = %d; ',                      data_FAMA.y_LEN);
fprintf(fw,'// size of the vector y\n');
fprintf(fw,'int LAM_LEN = %d; ',                    data_FAMA.lam_LEN);
fprintf(fw,'// size of the vector lam\n');
fprintf(fw,'int Q_LEN =  %d; ',                     data_FAMA.q_LEN);
fprintf(fw,'// size of the vector q\n');
fprintf(fw,'int T_LEN =  %d; ',                     data_FAMA.t_LEN);
fprintf(fw,'// size of the vector t\n');
fprintf(fw,'int CONST_VEC_LEN =  %d; ',             data_FAMA.const_vec_LEN);
fprintf(fw,'// size of the constant vector \n');
fprintf(fw,'int LX_LEN =  %d; ',                    data_FAMA.Lx_LEN);
fprintf(fw,'// storage for over relaxation \n');

% sizes of matrices
fprintf(fw,'int A_ROWS = %d; ',                     data_FAMA.A_ROWS);
fprintf(fw,'// number of rows in matrix A\n');
fprintf(fw,'int A_COLS = %d; ',                     data_FAMA.A_COLS);
fprintf(fw,'// number of columns in matrix A\n');
fprintf(fw,'int Q_ROWS = %d; ',                     data_FAMA.Q_ROWS);
fprintf(fw,'// number of rows in matrix Q\n');
fprintf(fw,'int Q_COLS = %d; ',                     data_FAMA.Q_COLS);
fprintf(fw,'// number of columns matrix Q\n');
fprintf(fw,'int L_ROWS = %d; ',                     data_FAMA.L_ROWS);
fprintf(fw,'// number of rows in matrix L\n');
fprintf(fw,'int L_COLS = %d; ',                     data_FAMA.L_COLS);
fprintf(fw,'// number of columns in matrix L\n');
fprintf(fw,'int SCALEDL_ROWS = %d; ',               data_FAMA.scaledL_ROWS);
fprintf(fw,'// number of rows in matrix rho*L\n');
fprintf(fw,'int SCALEDL_COLS = %d; ',               data_FAMA.scaledL_COLS);
fprintf(fw,'// number of cols in matrix rho*L\n');
fprintf(fw,'int MATRIX_ROWS =  %d; ',               data_FAMA.matrix_ROWS);
fprintf(fw,'// number of rows in matrix KKT\n');
fprintf(fw,'int MATRIX_COLS = %d; ',                data_FAMA.matrix_COLS);
fprintf(fw,'// number of cols in matrix KKT\n');
fprintf(fw,'int ALPHAL_ROWS =  %d; ',               data_FAMA.alphaL_ROWS);
fprintf(fw,'// number of rows in matrix alpha*L\n');
fprintf(fw,'int ALPHAL_COLS =  %d; ',               data_FAMA.alphaL_COLS);
fprintf(fw,'// number of cols in matrix alpha*L\n');


% sizes of matrices, containing indices of nonzero elements
fprintf(fw,'int  AI_ROWS =  %d; ',                  data_FAMA.AI_ROWS);
fprintf(fw,'// number of rows in matrix AI\n');
fprintf(fw,'int  AI_COLS =  %d; ',                  data_FAMA.AI_COLS);
fprintf(fw,'// number of columns in matrix AI\n');
fprintf(fw,'int  QI_ROWS =  %d; ',                  data_FAMA.QI_ROWS);
fprintf(fw,'// number of rows in matrix QI\n');
fprintf(fw,'int  QI_COLS =  %d; ',                  data_FAMA.QI_COLS);
fprintf(fw,'// number of columns matrix QI\n');
fprintf(fw,'int  LI_ROWS =  %d; ',                  data_FAMA.LI_ROWS);
fprintf(fw,'// number of rows in matrix LI\n');
fprintf(fw,'int  LI_COLS =  %d; ',                  data_FAMA.LI_COLS);
fprintf(fw,'// number of columns in matrix LI\n');
fprintf(fw,'int MATRIXI_ROWS =  %d; ',              data_FAMA.matrixI_ROWS);
fprintf(fw,'// number of rows in matrix KKT\n');
fprintf(fw,'int MATRIXI_COLS =  %d; ',              data_FAMA.matrixI_COLS);
fprintf(fw,'// number of cols in matrix KKT\n');
fprintf(fw,'int SCALEDLI_ROWS =  %d; ',             data_FAMA.scaledLI_ROWS);
fprintf(fw,'// number of rows in matrix rho*LI\n');
fprintf(fw,'int SCALEDLI_COLS =  %d; ',             data_FAMA.scaledLI_COLS);
fprintf(fw,'// number of cols in matrix rho*LI\n');
fprintf(fw,'int ALPHALI_ROWS =  %d; ',              data_FAMA.alphaLI_ROWS);
fprintf(fw,'// number of rows in matrix alpha*LI\n');
fprintf(fw,'int ALPHALI_COLS =  %d; ',              data_FAMA.alphaLI_COLS);
fprintf(fw,'// number of cols in matrix alpha*LI\n\n');

% export 1D arrays to a header file
generate_vector(filename, data_FAMA.x,              'x',              'double',  data_FAMA.x_LEN);
generate_vector(filename, data_FAMA.y,              'y',              'double',  data_FAMA.y_LEN);
generate_vector(filename, data_FAMA.prev_y,         'prev_y',         'double',  data_FAMA.y_LEN);
generate_vector(filename, data_FAMA.lam,            'lam',            'double',  data_FAMA.lam_LEN);
generate_vector(filename, data_FAMA.prev_lam,       'prev_lam',       'double',  data_FAMA.lam_LEN);
generate_vector(filename, data_FAMA.lam_hat,        'lam_hat',        'double',  data_FAMA.lam_LEN);
generate_vector(filename, data_FAMA.q,              'q',              'double',  data_FAMA.q_LEN);
generate_vector(filename, data_FAMA.b,              'b',              'double',  data_FAMA.b_LEN);
generate_vector(filename, data_FAMA.f,              'f',              'double',  data_FAMA.f_LEN);
generate_vector(filename, data_FAMA.c,              'c',              'double',  data_FAMA.c_LEN);
generate_vector(filename, data_FAMA.l,              'l',              'double',  data_FAMA.l_LEN);
generate_vector(filename, data_FAMA.const_vec,      'const_vec',      'double',  data_FAMA.const_vec_LEN);
generate_vector(filename, data_FAMA.new_prox_beg,   'new_prox_beg',   'int',     data_FAMA.nProx);
generate_vector(filename, data_FAMA.new_prox_end,   'new_prox_end',   'int',     data_FAMA.nProx);
generate_vector(filename, data_FAMA.t           ,   't',              'double',  data_FAMA.t_LEN);
generate_vector(filename, data_FAMA.lx,             'lx',             'double',  data_FAMA.Lx_LEN);
generate_vector(filename, data_FAMA.len_of_vectors, 'len_of_vectors', 'int',     data_FAMA.nProx);

% export matrices to a header file
generate_matrix(filename, data_FAMA.A,        'A',        'double',  data_FAMA.A_ROWS,        data_FAMA.A_COLS);
generate_matrix(filename, data_FAMA.Q,        'Q',        'double',  data_FAMA.Q_ROWS,        data_FAMA.Q_COLS);
generate_matrix(filename, data_FAMA.L,        'L',        'double',  data_FAMA.L_ROWS,        data_FAMA.L_COLS);
generate_matrix(filename, data_FAMA.matrix,   'matrix',   'double',  data_FAMA.matrix_ROWS,   data_FAMA.matrix_COLS);
generate_matrix(filename, data_FAMA.scaledL,  'scaledL',  'double',  data_FAMA.scaledL_ROWS,  data_FAMA.scaledL_COLS);
generate_matrix(filename, data_FAMA.alphaL,   'alphaL',   'double',  data_FAMA.alphaL_ROWS,   data_FAMA.alphaL_COLS);

% export matrice to a header file - indices of nonzero elements
generate_matrix(filename, data_FAMA.AI,        'AI',        'int',   data_FAMA.AI_ROWS,       data_FAMA.AI_COLS);
generate_matrix(filename, data_FAMA.QI,        'QI',        'int',   data_FAMA.QI_ROWS,       data_FAMA.QI_COLS);
generate_matrix(filename, data_FAMA.LI,        'LI',        'int',   data_FAMA.LI_ROWS,       data_FAMA.LI_COLS);
generate_matrix(filename, data_FAMA.matrixI,   'matrixI',   'int',   data_FAMA.matrixI_ROWS,  data_FAMA.matrixI_COLS);
generate_matrix(filename, data_FAMA.scaledLI,  'scaledLI',  'int',   data_FAMA.scaledLI_ROWS, data_FAMA.scaledLI_COLS);
generate_matrix(filename, data_FAMA.alphaLI,   'alphaLI',   'int',   data_FAMA.alphaLI_ROWS,  data_FAMA.alphaLI_COLS);

% generate array of proxes and norms
generate_array_of_prox_types(filename, data_FAMA.prox_type, data_FAMA.nProx);
generate_array_of_norm_types(filename, data_FAMA.normType,  data_FAMA.nProx);

% generate 1D arrays containing normTypesm, weights and constants to perform projections
generate_vector(filename, data_FAMA.weight,    'weights',   'double',  data_FAMA.nProx);
generate_vector(filename, data_FAMA.constants, 'constants', 'double',  data_FAMA.nProx);

% preprocessor directive ends
fw = fopen(filename, 'a');
fprintf(fw, '\n');
fprintf(fw, '#endif');

% close the file handler, otherwise std::exception is thrown
fclose(fw);
