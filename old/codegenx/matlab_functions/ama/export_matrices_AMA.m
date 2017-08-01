function export_matrices_AMA(data_AMA)

% create a file handler and the corresponding header file
filename = 'AMA_matrices.h';
fw = fopen(filename, 'w');

% do not allow mulitple includes
fprintf(fw, '#ifndef _AMA_MATRICES__H__\n');
fprintf(fw, '#define _AMA_MATRICES__H__\n\n');

% basic optimization data
fprintf(fw,'int NUM_OPT_VAR = %d; ',                data_AMA.num_of_opt_var);
fprintf(fw,'// number of optimisation variables\n');
fprintf(fw,'int NUM_OF_STATES = %d; ',              data_AMA.num_of_states);
fprintf(fw,'// number of states\n');
fprintf(fw,'int NPROX  = %d; ',                     data_AMA.nProx);
fprintf(fw,'// number of proxes\n');
fprintf(fw,'double DUALTOL =  %.15f; ',             data_AMA.settings.dualTol);
fprintf(fw,'// dual residual error\n');
fprintf(fw,'double PRIMALTOL = %.15f; ',            data_AMA.settings.primalTol);
fprintf(fw,'// primal residual error\n');
fprintf(fw,'double ALPHA = %.15f; ',                data_AMA.settings.relaxation);
fprintf(fw,'// parameter for over relaxation\n');
fprintf(fw,'int MAXITER  =  %d; ',                  data_AMA.settings.maxItr);
fprintf(fw,'// maximum number of iterations\n');
fprintf(fw,'int TERM_CHECK_FREQ =  %d; ',           data_AMA.settings.terminationCheckFreq);
fprintf(fw,'// frequency of termination check\n');
fprintf(fw,'double STEPSIZE = %.15f; ',             data_AMA.settings.rho);
fprintf(fw,'// stepsize\n');

% length of vectors
fprintf(fw,'int B_LEN  =  %d; ',                    data_AMA.b_LEN);
fprintf(fw,'// size of the vector b\n');
fprintf(fw,'int F_LEN =   %d; ',                    data_AMA.f_LEN);
fprintf(fw,'// size of the vector f\n'); 
fprintf(fw,'int C_LEN =  %d; ',                     data_AMA.c_LEN);
fprintf(fw,'// size of the vector c\n');
fprintf(fw,'int L_LEN = %d; ',                      data_AMA.l_LEN);
fprintf(fw,'// size of the vector l\n');
fprintf(fw,'int X_LEN = %d; ',                      data_AMA.x_LEN);
fprintf(fw,'// size of the vector x\n');
fprintf(fw,'int Y_LEN = %d; ',                      data_AMA.y_LEN);
fprintf(fw,'// size of the vector y\n');
fprintf(fw,'int LAM_LEN = %d; ',                    data_AMA.lam_LEN);
fprintf(fw,'// size of the vector lam\n');
fprintf(fw,'int Q_LEN =  %d; ',                     data_AMA.q_LEN);
fprintf(fw,'// size of the vector q\n');
fprintf(fw,'int T_LEN =  %d; ',                     data_AMA.t_LEN);
fprintf(fw,'// size of the vector t\n');
fprintf(fw,'int CONST_VEC_LEN =  %d; ',             data_AMA.const_vec_LEN);
fprintf(fw,'// size of the constant vector \n');
fprintf(fw,'int LX_LEN =  %d; ',                    data_AMA.Lx_LEN);
fprintf(fw,'// storage for over relaxation \n');

% sizes of matrices
fprintf(fw,'int A_ROWS = %d; ',                     data_AMA.A_ROWS);
fprintf(fw,'// number of rows in matrix A\n');
fprintf(fw,'int A_COLS = %d; ',                     data_AMA.A_COLS);
fprintf(fw,'// number of columns in matrix A\n');
fprintf(fw,'int Q_ROWS = %d; ',                     data_AMA.Q_ROWS);
fprintf(fw,'// number of rows in matrix Q\n');
fprintf(fw,'int Q_COLS = %d; ',                     data_AMA.Q_COLS);
fprintf(fw,'// number of columns matrix Q\n');
fprintf(fw,'int L_ROWS = %d; ',                     data_AMA.L_ROWS);
fprintf(fw,'// number of rows in matrix L\n');
fprintf(fw,'int L_COLS = %d; ',                     data_AMA.L_COLS);
fprintf(fw,'// number of columns in matrix L\n');
fprintf(fw,'int SCALEDL_ROWS = %d; ',               data_AMA.scaledL_ROWS);
fprintf(fw,'// number of rows in matrix rho*L\n');
fprintf(fw,'int SCALEDL_COLS = %d; ',               data_AMA.scaledL_COLS);
fprintf(fw,'// number of cols in matrix rho*L\n');
fprintf(fw,'int MATRIX_ROWS =  %d; ',               data_AMA.matrix_ROWS);
fprintf(fw,'// number of rows in matrix KKT\n');
fprintf(fw,'int MATRIX_COLS = %d; ',                data_AMA.matrix_COLS);
fprintf(fw,'// number of cols in matrix KKT\n');
fprintf(fw,'int ALPHAL_ROWS =  %d; ',               data_AMA.alphaL_ROWS);
fprintf(fw,'// number of rows in matrix alpha*LI\n');
fprintf(fw,'int ALPHAL_COLS =  %d; ',               data_AMA.alphaL_COLS);
fprintf(fw,'// number of cols in matrix alpha*LI\n');


% sizes of matrices, containing indices of nonzero elements
fprintf(fw,'int  AI_ROWS =  %d; ',                  data_AMA.AI_ROWS);
fprintf(fw,'// number of rows in matrix AI\n');
fprintf(fw,'int  AI_COLS =  %d; ',                  data_AMA.AI_COLS);
fprintf(fw,'// number of columns in matrix AI\n');
fprintf(fw,'int  QI_ROWS =  %d; ',                  data_AMA.QI_ROWS);
fprintf(fw,'// number of rows in matrix QI\n');
fprintf(fw,'int  QI_COLS =  %d; ',                  data_AMA.QI_COLS);
fprintf(fw,'// number of columns matrix QI\n');
fprintf(fw,'int  LI_ROWS =  %d; ',                  data_AMA.LI_ROWS);
fprintf(fw,'// number of rows in matrix LI\n');
fprintf(fw,'int  LI_COLS =  %d; ',                  data_AMA.LI_COLS);
fprintf(fw,'// number of columns in matrix LI\n');
fprintf(fw,'int MATRIXI_ROWS =  %d; ',              data_AMA.matrixI_ROWS);
fprintf(fw,'// number of rows in matrix KKT\n');
fprintf(fw,'int MATRIXI_COLS =  %d; ',              data_AMA.matrixI_COLS);
fprintf(fw,'// number of cols in matrix KKT\n');
fprintf(fw,'int SCALEDLI_ROWS =  %d; ',             data_AMA.scaledLI_ROWS);
fprintf(fw,'// number of rows in matrix rho*LI\n');
fprintf(fw,'int SCALEDLI_COLS =  %d; ',             data_AMA.scaledLI_COLS);
fprintf(fw,'// number of cols in matrix rho*LI\n');
fprintf(fw,'int ALPHALI_ROWS =  %d; ',              data_AMA.alphaLI_ROWS);
fprintf(fw,'// number of rows in matrix alpha*LI\n');
fprintf(fw,'int ALPHALI_COLS =  %d; ',              data_AMA.alphaLI_COLS);
fprintf(fw,'// number of cols in matrix alpha*LI\n\n');

% export 1D arrays to a header file
generate_vector(filename, data_AMA.x,              'x',              'double',  data_AMA.x_LEN);
generate_vector(filename, data_AMA.y,              'y',              'double',  data_AMA.y_LEN);
generate_vector(filename, data_AMA.prev_y,         'prev_y',         'double',  data_AMA.y_LEN);
generate_vector(filename, data_AMA.lam,            'lam',            'double',  data_AMA.lam_LEN);
generate_vector(filename, data_AMA.q,              'q',              'double',  data_AMA.q_LEN);
generate_vector(filename, data_AMA.b,              'b',              'double',  data_AMA.b_LEN);
generate_vector(filename, data_AMA.f,              'f',              'double',  data_AMA.f_LEN);
generate_vector(filename, data_AMA.c,              'c',              'double',  data_AMA.c_LEN);
generate_vector(filename, data_AMA.l,              'l',              'double',  data_AMA.l_LEN);
generate_vector(filename, data_AMA.const_vec,      'const_vec',      'double',  data_AMA.const_vec_LEN);
generate_vector(filename, data_AMA.new_prox_beg,   'new_prox_beg',   'int',     data_AMA.nProx);
generate_vector(filename, data_AMA.new_prox_end,   'new_prox_end',   'int',     data_AMA.nProx);
generate_vector(filename, data_AMA.t,              't',              'double',  data_AMA.t_LEN);
generate_vector(filename, data_AMA.lx,             'lx',             'double',  data_AMA.Lx_LEN);
generate_vector(filename, data_AMA.len_of_vectors, 'len_of_vectors', 'int',     data_AMA.nProx);

% export matrices to a header file
generate_matrix(filename, data_AMA.A,        'A',        'double',  data_AMA.A_ROWS,       data_AMA.A_COLS);
generate_matrix(filename, data_AMA.Q,        'Q',        'double',  data_AMA.Q_ROWS,       data_AMA.Q_COLS);
generate_matrix(filename, data_AMA.L,        'L',        'double',  data_AMA.L_ROWS,       data_AMA.L_COLS);
generate_matrix(filename, data_AMA.matrix,   'matrix',   'double',  data_AMA.matrix_ROWS,  data_AMA.matrix_COLS);
generate_matrix(filename, data_AMA.scaledL,  'scaledL',  'double',  data_AMA.scaledL_ROWS, data_AMA.scaledL_COLS);
generate_matrix(filename, data_AMA.alphaL,   'alphaL',   'double',  data_AMA.alphaL_ROWS,  data_AMA.alphaL_COLS);

% export matrice to a header file - indices of nonzero elements
generate_matrix(filename, data_AMA.AI,        'AI',        'int',  data_AMA.AI_ROWS,        data_AMA.AI_COLS);
generate_matrix(filename, data_AMA.QI,        'QI',        'int',  data_AMA.QI_ROWS,        data_AMA.QI_COLS);
generate_matrix(filename, data_AMA.LI,        'LI',        'int',  data_AMA.LI_ROWS,        data_AMA.LI_COLS);
generate_matrix(filename, data_AMA.matrixI,   'matrixI',   'int',  data_AMA.matrixI_ROWS,   data_AMA.matrixI_COLS);
generate_matrix(filename, data_AMA.scaledLI,  'scaledLI',  'int',  data_AMA.scaledLI_ROWS,  data_AMA.scaledLI_COLS);
generate_matrix(filename, data_AMA.alphaLI,   'alphaLI',   'int',  data_AMA.alphaLI_ROWS,   data_AMA.alphaLI_COLS);

% generate array of proxes and norms
generate_array_of_prox_types(filename, data_AMA.prox_type, data_AMA.nProx);
generate_array_of_norm_types(filename, data_AMA.normType,  data_AMA.nProx);

% generate 1D arrays containing normTypesm, weights and constants to perform projections
generate_vector(filename, data_AMA.weight,    'weights',   'double',  data_AMA.nProx);
generate_vector(filename, data_AMA.constants, 'constants', 'double',  data_AMA.nProx);

% preprocessor directive ends
fw = fopen(filename, 'a');
fprintf(fw, '\n');
fprintf(fw, '#endif');

% close the file handler, otherwise std::exception is thrown
fclose(fw);