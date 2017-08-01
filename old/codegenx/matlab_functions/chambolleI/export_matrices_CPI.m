function export_matrices_CPI(data_CPI)

% create a file handler and the corresponding header file
filename = 'CPI_matrices.h';
fw = fopen(filename, 'w');

% do not allow mulitple includes
fprintf(fw, '#ifndef _CPI_MATRICES__H__\n');
fprintf(fw, '#define _CPI_MATRICES__H__\n\n');

% basic optimization data
fprintf(fw,'int NUM_OPT_VAR = %d; ',                data_CPI.num_of_opt_var);
fprintf(fw,'// number of optimisation variables\n');
fprintf(fw,'int NUM_OF_STATES = %d; ',              data_CPI.num_of_states);
fprintf(fw,'// number of states\n');
fprintf(fw,'int NPROX  = %d; ',                     data_CPI.nProx);
fprintf(fw,'// number of proxes\n');
fprintf(fw,'double DUALTOL =  %.15f; ',             data_CPI.settings.dualTol);
fprintf(fw,'// dual residual error\n');
fprintf(fw,'double PRIMALTOL = %.15f; ',            data_CPI.settings.primalTol);
fprintf(fw,'// primal residual error\n');
fprintf(fw,'int MAXITER  =  %d; ',                  data_CPI.settings.maxItr);
fprintf(fw,'// maximum number of iterations\n');
fprintf(fw,'int TERM_CHECK_FREQ =  %d; ',           data_CPI.settings.terminationCheckFreq);
fprintf(fw,'//  frequency of termination check\n');
fprintf(fw,'double THETA  =  %.15f; ',              data_CPI.theta);
fprintf(fw,'// maximum number of iterations\n');
fprintf(fw,'double SIGMA  =  %.15f; ',              data_CPI.sigma);
fprintf(fw,'// maximum number of iterations\n');
fprintf(fw,'double TAU  =  %.15f; ',                data_CPI.tau);
fprintf(fw,'// maximum number of iterations\n');

% length of vectors
fprintf(fw,'int B_LEN  =  %d; ',                    data_CPI.b_LEN);
fprintf(fw,'// size of the vector b\n');
fprintf(fw,'int F_LEN =   %d; ',                    data_CPI.f_LEN);
fprintf(fw,'// size of the vector f\n'); 
fprintf(fw,'int C_LEN =  %d; ',                     data_CPI.c_LEN);
fprintf(fw,'// size of the vector c\n');
fprintf(fw,'int L_LEN = %d; ',                      data_CPI.l_LEN);
fprintf(fw,'// size of the vector l\n');
fprintf(fw,'int X_LEN = %d; ',                      data_CPI.x_LEN);
fprintf(fw,'// size of the vector x\n');
fprintf(fw,'int NU_LEN = %d; ',                     data_CPI.nu_LEN);
fprintf(fw,'// size of the vector nu\n');
fprintf(fw,'int P_LEN = %d; ',                      data_CPI.p_LEN);
fprintf(fw,'// size of the vector p\n');
fprintf(fw,'int P_NU_LEN =  %d; ',                  data_CPI.p_nu_LEN);
fprintf(fw,'// size of the vector p_nu\n');

% sizes of matrices
fprintf(fw,'int A_ROWS = %d; ',                     data_CPI.A_ROWS);
fprintf(fw,'// number of rows in matrix A\n');
fprintf(fw,'int A_COLS = %d; ',                     data_CPI.A_COLS);
fprintf(fw,'// number of columns in matrix A\n');
fprintf(fw,'int Q_ROWS = %d; ',                     data_CPI.Q_ROWS);
fprintf(fw,'// number of rows in matrix Q\n');
fprintf(fw,'int Q_COLS = %d; ',                     data_CPI.Q_COLS);
fprintf(fw,'// number of columns matrix Q\n');
fprintf(fw,'int L_ROWS = %d; ',                     data_CPI.L_ROWS);
fprintf(fw,'// number of rows in matrix L\n');
fprintf(fw,'int L_COLS = %d; ',                     data_CPI.L_COLS);
fprintf(fw,'// number of columns in matrix L\n');
fprintf(fw,'int K_ROWS = %d; ',                     data_CPI.K_ROWS);
fprintf(fw,'// number of rows in matrix K\n');
fprintf(fw,'int K_COLS = %d; ',                     data_CPI.K_COLS);
fprintf(fw,'// number of columns in matrix K\n');
fprintf(fw,'int KT_ROWS = %d; ',                    data_CPI.KT_ROWS);
fprintf(fw,'// number of rows in matrix K\n');
fprintf(fw,'int KT_COLS = %d; ',                    data_CPI.KT_COLS);
fprintf(fw,'// number of columns in matrix K\n');
fprintf(fw,'int MINVERSE_ROWS =  %d; ',             data_CPI.mInverse_ROWS);
fprintf(fw,'// number of rows in inverse_matrix\n');
fprintf(fw,'int MINVERSE_COLS = %d; ',              data_CPI.mInverse_COLS);
fprintf(fw,'// number of cols in in inverse_matrix\n');
fprintf(fw,'int SP_ROWS = %d; ',                    data_CPI.SP_ROWS);
fprintf(fw,'// number of rows in matrix SP\n');
fprintf(fw,'int SP_COLS = %d; ',                    data_CPI.SP_COLS);
fprintf(fw,'// number of cols in matrix SP\n');
fprintf(fw,'int SD1_ROWS = %d; ',                   data_CPI.SD1_ROWS);
fprintf(fw,'// number of rows in inverse_matrix\n');
fprintf(fw,'int SD1_COLS = %d; ',                   data_CPI.SD1_COLS);
fprintf(fw,'// number of cols in in inverse_matrix\n');
fprintf(fw,'int SD2_ROWS = %d; ',                   data_CPI.SD2_ROWS);
fprintf(fw,'// number of rows in inverse_matrix\n');
fprintf(fw,'int SD2_COLS = %d; ',                   data_CPI.SD2_COLS);
fprintf(fw,'// number of cols in in inverse_matrix\n');

% sizes of matrices, containing indices of nonzero elements
fprintf(fw,'int  AI_ROWS =  %d; ',                  data_CPI.AI_ROWS);
fprintf(fw,'// number of rows in matrix AI\n');
fprintf(fw,'int  AI_COLS =  %d; ',                  data_CPI.AI_COLS);
fprintf(fw,'// number of columns in matrix AI\n');
fprintf(fw,'int  QI_ROWS =  %d; ',                  data_CPI.QI_ROWS);
fprintf(fw,'// number of rows in matrix QI\n');
fprintf(fw,'int  QI_COLS =  %d; ',                  data_CPI.QI_COLS);
fprintf(fw,'// number of columns matrix QI\n');
fprintf(fw,'int  KI_ROWS =  %d; ',                  data_CPI.KI_ROWS);
fprintf(fw,'// number of rows in matrix KI\n');
fprintf(fw,'int  KI_COLS =  %d; ',                  data_CPI.KI_COLS);
fprintf(fw,'// number of columns in matrix KI\n');
fprintf(fw,'int  KTI_ROWS =  %d; ',                 data_CPI.KTI_ROWS);
fprintf(fw,'// number of rows in matrix KI\n');
fprintf(fw,'int  KTI_COLS =  %d; ',                 data_CPI.KTI_COLS);
fprintf(fw,'// number of columns in matrix KI\n');
fprintf(fw,'int  LI_ROWS =  %d; ',                  data_CPI.LI_ROWS);
fprintf(fw,'// number of rows in matrix LI\n');
fprintf(fw,'int  LI_COLS =  %d; ',                  data_CPI.LI_COLS);
fprintf(fw,'// number of columns in matrix LI\n');
fprintf(fw,'int MINVERSEI_COLS =  %d; ',            data_CPI.mInverseI_ROWS);
fprintf(fw,'// number of rows in inverse_matrix\n');
fprintf(fw,'int MINVERSEI_ROWS =  %d; ',            data_CPI.mInverseI_COLS);
fprintf(fw,'// number of cols in inverse_matrix\n');
fprintf(fw,'int SPI_ROWS =  %d; ',                  data_CPI.SPI_ROWS);
fprintf(fw,'// number of rows in matrix SPI\n');
fprintf(fw,'int SPI_COLS =  %d; ',                  data_CPI.SPI_COLS);
fprintf(fw,'// number of cols in matrix SPI\n');
fprintf(fw,'int SD1I_ROWS = %d; ',                  data_CPI.SD1I_ROWS);
fprintf(fw,'// number of rows in matrix SD1I\n');
fprintf(fw,'int SD1I_COLS = %d; ',                  data_CPI.SD1I_COLS);
fprintf(fw,'// number of cols in matrix SD1I\n');
fprintf(fw,'int SD2I_ROWS = %d; ',                  data_CPI.SD2I_ROWS);
fprintf(fw,'// number of rows in matrix SD2I\n');
fprintf(fw,'int SD2I_COLS = %d; ',                  data_CPI.SD2I_COLS);
fprintf(fw,'// number of cols in matrix SD2I\n\n');

% export 1D arrays to a header file
generate_vector(filename, data_CPI.x,              'x',              'double',  data_CPI.x_LEN);
generate_vector(filename, data_CPI.x_prev,         'x_prev',         'double',  data_CPI.x_LEN);
generate_vector(filename, data_CPI.x,              'xbar',           'double',  data_CPI.x_LEN);
generate_vector(filename, data_CPI.xbar_prev,      'xbar_prev',      'double',  data_CPI.x_LEN);
generate_vector(filename, data_CPI.nu,             'nu',             'double',  data_CPI.nu_LEN);
generate_vector(filename, data_CPI.nu_prev,        'nu_prev',        'double',  data_CPI.nu_LEN);
generate_vector(filename, data_CPI.p,              'p',              'double',  data_CPI.p_LEN);
generate_vector(filename, data_CPI.p_prev,         'p_prev',         'double',  data_CPI.p_LEN);
generate_vector(filename, data_CPI.p_nu,           'p_nu',           'double',  data_CPI.p_nu_LEN);
generate_vector(filename, data_CPI.p_nu_prev,      'p_nu_prev',      'double',  data_CPI.p_nu_LEN);
generate_vector(filename, data_CPI.b,              'b',              'double',  data_CPI.b_LEN);
generate_vector(filename, data_CPI.f,              'f',              'double',  data_CPI.f_LEN);
generate_vector(filename, data_CPI.c,              'c',              'double',  data_CPI.c_LEN);
generate_vector(filename, data_CPI.l,              'l',              'double',  data_CPI.l_LEN);
generate_vector(filename, data_CPI.new_prox_beg,   'new_prox_beg',   'int',     data_CPI.nProx);
generate_vector(filename, data_CPI.new_prox_end,   'new_prox_end',   'int',     data_CPI.nProx);
generate_vector(filename, data_CPI.len_of_vectors, 'len_of_vectors', 'int',     data_CPI.nProx);

% export matrices to a header file
generate_matrix(filename, data_CPI.A,        'A',        'double',  data_CPI.A_ROWS,        data_CPI.A_COLS);
generate_matrix(filename, data_CPI.Q,        'Q',        'double',  data_CPI.Q_ROWS,        data_CPI.Q_COLS);
generate_matrix(filename, data_CPI.K,        'K',        'double',  data_CPI.K_ROWS,        data_CPI.K_COLS);
generate_matrix(filename, data_CPI.KT,       'KT',       'double',  data_CPI.KT_ROWS,       data_CPI.KT_COLS);
generate_matrix(filename, data_CPI.L,        'L',        'double',  data_CPI.L_ROWS,        data_CPI.L_COLS);
generate_matrix(filename, data_CPI.mInverse, 'mInverse', 'double',  data_CPI.mInverse_ROWS, data_CPI.mInverse_COLS);
generate_matrix(filename, data_CPI.SP,       'SP',       'double',  data_CPI.SP_ROWS,       data_CPI.SP_COLS);
generate_matrix(filename, data_CPI.SD1,      'SD1',      'double',  data_CPI.SD1_ROWS,      data_CPI.SD1_COLS);
generate_matrix(filename, data_CPI.SD2,      'SD2',      'double',  data_CPI.SD2_ROWS,      data_CPI.SD2_COLS);

% export matrice to a header file - indices of nonzero elements
generate_matrix(filename, data_CPI.AI,        'AI',        'int',  data_CPI.AI_ROWS,        data_CPI.AI_COLS);
generate_matrix(filename, data_CPI.QI,        'QI',        'int',  data_CPI.QI_ROWS,        data_CPI.QI_COLS);
generate_matrix(filename, data_CPI.KI,        'KI',        'int',  data_CPI.KI_ROWS,        data_CPI.KI_COLS);
generate_matrix(filename, data_CPI.KTI,       'KTI',       'int',  data_CPI.KTI_ROWS,       data_CPI.KTI_COLS);
generate_matrix(filename, data_CPI.LI,        'LI',        'int',  data_CPI.LI_ROWS,        data_CPI.LI_COLS);
generate_matrix(filename, data_CPI.mInverseI, 'mInverseI', 'int',  data_CPI.mInverseI_ROWS, data_CPI.mInverseI_COLS);
generate_matrix(filename, data_CPI.SPI,       'SPI',       'int',  data_CPI.SPI_ROWS,       data_CPI.SPI_COLS);
generate_matrix(filename, data_CPI.SD1I,      'SD1I',      'int',  data_CPI.SD1I_ROWS,      data_CPI.SD1I_COLS);
generate_matrix(filename, data_CPI.SD2I,      'SD2I',      'int',  data_CPI.SD2I_ROWS,      data_CPI.SD2I_COLS);

% generate array of proxes and norms
generate_array_of_prox_types(filename, data_CPI.proxConj_type, data_CPI.nProx);
generate_array_of_norm_types(filename, data_CPI.dualNormType,  data_CPI.nProx);

% generate 1D arrays containing normTypesm, weights and constants to perform projections
generate_vector(filename, data_CPI.weight,    'weights',   'double',  data_CPI.nProx);
generate_vector(filename, data_CPI.constants, 'constants', 'double',  data_CPI.nProx);

% preprocessor directive ends
fw = fopen(filename, 'a');
fprintf(fw, '\n');
fprintf(fw, '#endif');

% close the file handler, otherwise std::exception is thrown
fclose(fw);