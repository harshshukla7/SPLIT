function export_matrices_CPIprec(data_CPIprec)

% create a file handler and the corresponding header file
filename = 'CHAMBOLLEI_matrices.h';
fw = fopen(filename, 'w');

% do not allow mulitple includes
fprintf(fw, '#ifndef _CHAMBOLLEI_MATRICES__H__\n');
fprintf(fw, '#define _CHAMBOLLEI_MATRICES__H__\n\n');

% basic optimization data
fprintf(fw,'int NUM_OPT_VAR = %d; ',                data_CPIprec.num_of_opt_var);
fprintf(fw,'// number of optimisation variables\n');
fprintf(fw,'int NUM_OF_STATES = %d; ',              data_CPIprec.num_of_states);
fprintf(fw,'// number of states\n');
fprintf(fw,'int NPROX  = %d; ',                     data_CPIprec.nProx);
fprintf(fw,'// number of proxes\n');
fprintf(fw,'float DUALTOL =  %f; ',                 data_CPIprec.settings.dualTol);
fprintf(fw,'// dual residual error\n');
fprintf(fw,'float PRIMALTOL = %f; ',                data_CPIprec.settings.primalTol);
fprintf(fw,'// primal residual error\n');
fprintf(fw,'int MAXITER  =  %d; ',                  data_CPIprec.settings.maxItr);
fprintf(fw,'// maximum number of iterations\n');
fprintf(fw,'float STEPSIZE = %f; ',                 data_CPIprec.settings.rho);
fprintf(fw,'// value of the stepsize\n');
fprintf(fw,'int TERM_CHECK_FREQ =  %d; ',           data_CPIprec.settings.terminationCheckFreq);
fprintf(fw,'//  frequency of termination check\n');
fprintf(fw,'float THETA  =  %f; ',                  data_CPIprec.theta);
fprintf(fw,'// maximum number of iterations\n');
fprintf(fw,'float SIGMA  =  %f; ',                  data_CPIprec.sigma);
fprintf(fw,'// maximum number of iterations\n');
fprintf(fw,'float TAU  =  %f; ',                    data_CPIprec.tau);
fprintf(fw,'// maximum number of iterations\n');
fprintf(fw,'float TAUTHETA  =  %f; ',               data_CPIprec.invtauTheta);
fprintf(fw,'// maximum number of iterations\n');

% length of vectors
fprintf(fw,'int B_LEN  =  %d; ',                    data_CPIprec.b_LEN);
fprintf(fw,'// size of the vector b\n');
fprintf(fw,'int F_LEN =   %d; ',                    data_CPIprec.f_LEN);
fprintf(fw,'// size of the vector f\n'); 
fprintf(fw,'int C_LEN =  %d; ',                     data_CPIprec.c_LEN);
fprintf(fw,'// size of the vector c\n');
fprintf(fw,'int L_LEN = %d; ',                      data_CPIprec.l_LEN);
fprintf(fw,'// size of the vector l\n');
fprintf(fw,'int X_LEN = %d; ',                      data_CPIprec.x_LEN);
fprintf(fw,'// size of the vector x\n');
fprintf(fw,'int NU_LEN = %d; ',                     data_CPIprec.nu_LEN);
fprintf(fw,'// size of the vector nu\n');
fprintf(fw,'int P_LEN = %d; ',                      data_CPIprec.p_LEN);
fprintf(fw,'// size of the vector p\n');
fprintf(fw,'int P_NU_LEN =  %d; ',                  data_CPIprec.p_nu_LEN);
fprintf(fw,'// size of the vector p_nu\n');

% sizes of matrices
fprintf(fw,'int A_ROWS = %d; ',                     data_CPIprec.A_ROWS);
fprintf(fw,'// number of rows in matrix A\n');
fprintf(fw,'int A_COLS = %d; ',                     data_CPIprec.A_COLS);
fprintf(fw,'// number of columns in matrix A\n');
fprintf(fw,'int Q_ROWS = %d; ',                     data_CPIprec.Q_ROWS);
fprintf(fw,'// number of rows in matrix Q\n');
fprintf(fw,'int Q_COLS = %d; ',                     data_CPIprec.Q_COLS);
fprintf(fw,'// number of columns matrix Q\n');
fprintf(fw,'int L_ROWS = %d; ',                     data_CPIprec.L_ROWS);
fprintf(fw,'// number of rows in matrix L\n');
fprintf(fw,'int L_COLS = %d; ',                     data_CPIprec.L_COLS);
fprintf(fw,'// number of columns in matrix L\n');
fprintf(fw,'int K_ROWS = %d; ',                     data_CPIprec.K_ROWS);
fprintf(fw,'// number of rows in matrix K\n');
fprintf(fw,'int K_COLS = %d; ',                     data_CPIprec.K_COLS);
fprintf(fw,'// number of columns in matrix K\n');
fprintf(fw,'int KT_ROWS = %d; ',                    data_CPIprec.KT_ROWS);
fprintf(fw,'// number of rows in matrix K\n');
fprintf(fw,'int KT_COLS = %d; ',                    data_CPIprec.KT_COLS);
fprintf(fw,'// number of columns in matrix K\n');
fprintf(fw,'int SCALEDK_ROWS = %d; ',               data_CPIprec.scaledK_ROWS);
fprintf(fw,'// number of rows in matrix theta*K:\n');
fprintf(fw,'int SCALEDK_COLS = %d; ',               data_CPIprec.scaledK_COLS);
fprintf(fw,'// number of cols in matrix theta*K\n');
fprintf(fw,'int MINVERSE_ROWS =  %d; ',             data_CPIprec.mInverse_ROWS);
fprintf(fw,'// number of rows in inverse_matrix\n');
fprintf(fw,'int MINVERSE_COLS = %d; ',              data_CPIprec.mInverse_COLS);
fprintf(fw,'// number of cols in in inverse_matrix\n');

% sizes of matrices, containing indices of nonzero elements
fprintf(fw,'int  AI_ROWS =  %d; ',                  data_CPIprec.AI_ROWS);
fprintf(fw,'// number of rows in matrix AI\n');
fprintf(fw,'int  AI_COLS =  %d; ',                  data_CPIprec.AI_COLS);
fprintf(fw,'// number of columns in matrix AI\n');
fprintf(fw,'int  QI_ROWS =  %d; ',                  data_CPIprec.QI_ROWS);
fprintf(fw,'// number of rows in matrix QI\n');
fprintf(fw,'int  QI_COLS =  %d; ',                  data_CPIprec.QI_COLS);
fprintf(fw,'// number of columns matrix QI\n');
fprintf(fw,'int  KI_ROWS =  %d; ',                  data_CPIprec.KI_ROWS);
fprintf(fw,'// number of rows in matrix KI\n');
fprintf(fw,'int  KI_COLS =  %d; ',                  data_CPIprec.KI_COLS);
fprintf(fw,'// number of columns in matrix KI\n');
fprintf(fw,'int  KTI_ROWS =  %d; ',                 data_CPIprec.KTI_ROWS);
fprintf(fw,'// number of rows in matrix KI\n');
fprintf(fw,'int  KTI_COLS =  %d; ',                 data_CPIprec.KTI_COLS);
fprintf(fw,'// number of columns in matrix KI\n');
fprintf(fw,'int  LI_ROWS =  %d; ',                  data_CPIprec.LI_ROWS);
fprintf(fw,'// number of rows in matrix LI\n');
fprintf(fw,'int  LI_COLS =  %d; ',                  data_CPIprec.LI_COLS);
fprintf(fw,'// number of columns in matrix LI\n');
fprintf(fw,'int MINVERSEI_COLS =  %d; ',            data_CPIprec.mInverseI_ROWS);
fprintf(fw,'// number of rows in inverse_matrix\n');
fprintf(fw,'int MINVERSEI_ROWS =  %d; ',            data_CPIprec.mInverseI_COLS);
fprintf(fw,'// number of cols in inverse_matrix\n');
fprintf(fw,'int SCALEDKI_ROWS =  %d; ',             data_CPIprec.scaledKI_ROWS);
fprintf(fw,'// number of rows in matrix theta*K\n');
fprintf(fw,'int SCALEDKI_COLS =  %d; ',             data_CPIprec.scaledKI_COLS);
fprintf(fw,'// number of cols in matrix theta*K\n\n');

% export 1D arrays to a header file
generate_vector(filename, data_CPIprec.x,              'x',              'float',  data_CPIprec.x_LEN);
generate_vector(filename, data_CPIprec.x_prev,         'x_prev',         'float',  data_CPIprec.x_LEN);
generate_vector(filename, data_CPIprec.x,              'xbar',           'float',  data_CPIprec.x_LEN);
generate_vector(filename, data_CPIprec.nu,             'nu',             'float',  data_CPIprec.nu_LEN);
generate_vector(filename, data_CPIprec.nu_prev,        'nu_prev',        'float',  data_CPIprec.nu_LEN);
generate_vector(filename, data_CPIprec.p,              'p',              'float',  data_CPIprec.p_LEN);
generate_vector(filename, data_CPIprec.p_prev,         'p_prev',         'float',  data_CPIprec.p_LEN);
generate_vector(filename, data_CPIprec.p_nu,           'p_nu',           'float',  data_CPIprec.p_nu_LEN);
generate_vector(filename, data_CPIprec.p_nu_prev,      'p_nu_prev',      'float',  data_CPIprec.p_nu_LEN);
generate_vector(filename, data_CPIprec.b,              'b',              'float',  data_CPIprec.b_LEN);
generate_vector(filename, data_CPIprec.f,              'f',              'float',  data_CPIprec.f_LEN);
generate_vector(filename, data_CPIprec.c,              'c',              'float',  data_CPIprec.c_LEN);
generate_vector(filename, data_CPIprec.l,              'l',              'float',  data_CPIprec.l_LEN);
generate_vector(filename, data_CPIprec.new_prox_beg,   'new_prox_beg',   'int',    data_CPIprec.nProx);
generate_vector(filename, data_CPIprec.new_prox_end,   'new_prox_end',   'int',    data_CPIprec.nProx);
generate_vector(filename, data_CPIprec.len_of_vectors, 'len_of_vectors', 'int',    data_CPIprec.nProx);

% export matrices to a header file
generate_matrix(filename, data_CPIprec.A,        'A',        'float',  data_CPIprec.A_ROWS,        data_CPIprec.A_COLS);
generate_matrix(filename, data_CPIprec.Q,        'Q',        'float',  data_CPIprec.Q_ROWS,        data_CPIprec.Q_COLS);
generate_matrix(filename, data_CPIprec.K,        'K',        'float',  data_CPIprec.K_ROWS,        data_CPIprec.K_COLS);
generate_matrix(filename, data_CPIprec.KT,       'KT',       'float',  data_CPIprec.KT_ROWS,       data_CPIprec.KT_COLS);
generate_matrix(filename, data_CPIprec.L,        'L',        'float',  data_CPIprec.L_ROWS,        data_CPIprec.L_COLS);
generate_matrix(filename, data_CPIprec.mInverse, 'mInverse', 'float',  data_CPIprec.mInverse_ROWS, data_CPIprec.mInverse_COLS);
generate_matrix(filename, data_CPIprec.scaledK,  'scaledK',  'float',  data_CPIprec.scaledK_ROWS,  data_CPIprec.scaledK_COLS);

% export matrice to a header file - indices of nonzero elements
generate_matrix(filename, data_CPIprec.AI,        'AI',        'int',  data_CPIprec.AI_ROWS,        data_CPIprec.AI_COLS);
generate_matrix(filename, data_CPIprec.QI,        'QI',        'int',  data_CPIprec.QI_ROWS,        data_CPIprec.QI_COLS);
generate_matrix(filename, data_CPIprec.KI,        'KI',        'int',  data_CPIprec.KI_ROWS,        data_CPIprec.KI_COLS);
generate_matrix(filename, data_CPIprec.KTI,       'KTI',       'int',  data_CPIprec.KTI_ROWS,       data_CPIprec.KTI_COLS);
generate_matrix(filename, data_CPIprec.LI,        'LI',        'int',  data_CPIprec.LI_ROWS,        data_CPIprec.LI_COLS);
generate_matrix(filename, data_CPIprec.mInverseI, 'mInverseI', 'int',  data_CPIprec.mInverseI_ROWS, data_CPIprec.mInverseI_COLS);
generate_matrix(filename, data_CPIprec.scaledKI,  'scaledKI',  'int',  data_CPIprec.scaledKI_ROWS,  data_CPIprec.scaledKI_COLS);

% generate array of proxes and norms
generate_array_of_prox_types(filename, data_CPIprec.proxConj_type, data_CPIprec.nProx);
generate_array_of_norm_types(filename, data_CPIprec.dualNormType,  data_CPIprec.nProx);

% generate 1D arrays containing normTypesm, weights and constants to perform projections
generate_vector(filename, data_CPIprec.weight,    'weights',   'float',  data_CPIprec.nProx);
generate_vector(filename, data_CPIprec.constants, 'constants', 'float',  data_CPIprec.nProx);

% preprocessor directive ends
fw = fopen(filename, 'a');
fprintf(fw, '\n');
fprintf(fw, '#endif');

% close the file handler, otherwise std::exception is thrown
fclose(fw);