function export_matrices_CPII(data_CPII)

% create a file handler and the corresponding header file
filename = 'CPII_matrices.h';
fw = fopen(filename, 'w');

% do not allow mulitple includes
fprintf(fw, '#ifndef _CPII_MATRICES__H__\n');
fprintf(fw, '#define _CPII_MATRICES__H__\n\n');

% basic optimization data
fprintf(fw,'int NUM_OPT_VAR = %d; ',                data_CPII.num_of_opt_var);
fprintf(fw,'// number of optimisation variables\n');
fprintf(fw,'int NUM_OF_STATES = %d; ',              data_CPII.num_of_states);
fprintf(fw,'// number of states\n');
fprintf(fw,'int NPROX  = %d; ',                     data_CPII.nProx);
fprintf(fw,'// number of proxes\n');
fprintf(fw,'double DUALTOL =  %.15f; ',             data_CPII.settings.dualTol);
fprintf(fw,'// dual residual error\n');
fprintf(fw,'double PRIMALTOL = %.15f; ',            data_CPII.settings.primalTol);
fprintf(fw,'// primal residual error\n');
fprintf(fw,'int MAXITER  =  %d; ',                  data_CPII.settings.maxItr);
fprintf(fw,'// maximum number of iterations\n');
fprintf(fw,'int TERM_CHECK_FREQ =  %d; ',           data_CPII.settings.terminationCheckFreq);
fprintf(fw,'//  frequency of termination check\n');
fprintf(fw,'double THETA  =  %.15f; ',              data_CPII.theta);
fprintf(fw,'// stepsize in chambolle \n');
fprintf(fw,'double SIGMA  =  %.15f; ',              data_CPII.sigma);
fprintf(fw,'// stepsize in chambolle \n');
fprintf(fw,'double TAU  =  %.15f; ',                data_CPII.tau);
fprintf(fw,'// stepsize in chambolle \n');
fprintf(fw,'double GAMMA  =  %.15f; ',              data_CPII.gamma);
fprintf(fw,'// stepsize in chambolle \n');

% length of vectors
fprintf(fw,'int B_LEN  =  %d; ',                    data_CPII.b_LEN);
fprintf(fw,'// size of the vector b\n');
fprintf(fw,'int F_LEN =   %d; ',                    data_CPII.f_LEN);
fprintf(fw,'// size of the vector f\n'); 
fprintf(fw,'int C_LEN =  %d; ',                     data_CPII.c_LEN);
fprintf(fw,'// size of the vector c\n');
fprintf(fw,'int L_LEN = %d; ',                      data_CPII.l_LEN);
fprintf(fw,'// size of the vector l\n');
fprintf(fw,'int X_LEN = %d; ',                      data_CPII.x_LEN);
fprintf(fw,'// size of the vector x\n');
fprintf(fw,'int NU_LEN = %d; ',                     data_CPII.nu_LEN);
fprintf(fw,'// size of the vector nu\n');
fprintf(fw,'int P_LEN = %d; ',                      data_CPII.p_LEN);
fprintf(fw,'// size of the vector p\n');
fprintf(fw,'int P_NU_LEN =  %d; ',                  data_CPII.p_nu_LEN);
fprintf(fw,'// size of the vector p_nu\n');
fprintf(fw,'int DIAGQ_LEN = %d; ',                  data_CPII.diagQ_LEN);
fprintf(fw,'// size of the vector diagQ\n');

% sizes of matrices
fprintf(fw,'int A_ROWS = %d; ',                     data_CPII.A_ROWS);
fprintf(fw,'// number of rows in matrix A\n');
fprintf(fw,'int A_COLS = %d; ',                     data_CPII.A_COLS);
fprintf(fw,'// number of columns in matrix A\n');
fprintf(fw,'int L_ROWS = %d; ',                     data_CPII.L_ROWS);
fprintf(fw,'// number of rows in matrix L\n');
fprintf(fw,'int L_COLS = %d; ',                     data_CPII.L_COLS);
fprintf(fw,'// number of columns in matrix L\n');
fprintf(fw,'int K_ROWS = %d; ',                     data_CPII.K_ROWS);
fprintf(fw,'// number of rows in matrix K\n');
fprintf(fw,'int K_COLS = %d; ',                     data_CPII.K_COLS);
fprintf(fw,'// number of columns in matrix K\n');
fprintf(fw,'int KT_ROWS = %d; ',                    data_CPII.KT_ROWS);
fprintf(fw,'// number of rows in matrix K\n');
fprintf(fw,'int KT_COLS = %d; ',                    data_CPII.KT_COLS);
fprintf(fw,'// number of columns in matrix K\n');
fprintf(fw,'int SP_ROWS = %d; ',                    data_CPII.SP_ROWS);
fprintf(fw,'// number of rows in matrix SP\n');
fprintf(fw,'int SP_COLS = %d; ',                    data_CPII.SP_COLS);
fprintf(fw,'// number of cols in matrix SP\n');
fprintf(fw,'int SD1_ROWS = %d; ',                   data_CPII.SD1_ROWS);
fprintf(fw,'// number of rows in inverse_matrix\n');
fprintf(fw,'int SD1_COLS = %d; ',                   data_CPII.SD1_COLS);
fprintf(fw,'// number of cols in in inverse_matrix\n');
fprintf(fw,'int SD2_ROWS = %d; ',                   data_CPII.SD2_ROWS);
fprintf(fw,'// number of rows in inverse_matrix\n');
fprintf(fw,'int SD2_COLS = %d; ',                   data_CPII.SD2_COLS);
fprintf(fw,'// number of cols in in inverse_matrix\n');

% sizes of matrices, containing indices of nonzero elements
fprintf(fw,'int  AI_ROWS =  %d; ',                  data_CPII.AI_ROWS);
fprintf(fw,'// number of rows in matrix AI\n');
fprintf(fw,'int  AI_COLS =  %d; ',                  data_CPII.AI_COLS);
fprintf(fw,'// number of columns in matrix AI\n');
fprintf(fw,'int  KI_ROWS =  %d; ',                  data_CPII.KI_ROWS);
fprintf(fw,'// number of rows in matrix KI\n');
fprintf(fw,'int  KI_COLS =  %d; ',                  data_CPII.KI_COLS);
fprintf(fw,'// number of columns in matrix KI\n');
fprintf(fw,'int  KTI_ROWS =  %d; ',                 data_CPII.KTI_ROWS);
fprintf(fw,'// number of rows in matrix KI\n');
fprintf(fw,'int  KTI_COLS =  %d; ',                 data_CPII.KTI_COLS);
fprintf(fw,'// number of columns in matrix KI\n');
fprintf(fw,'int  LI_ROWS =  %d; ',                  data_CPII.LI_ROWS);
fprintf(fw,'// number of rows in matrix LI\n');
fprintf(fw,'int  LI_COLS =  %d; ',                  data_CPII.LI_COLS);
fprintf(fw,'// number of columns in matrix LI\n');
fprintf(fw,'int SPI_ROWS =  %d; ',                  data_CPII.SPI_ROWS);
fprintf(fw,'// number of rows in matrix SPI\n');
fprintf(fw,'int SPI_COLS =  %d; ',                  data_CPII.SPI_COLS);
fprintf(fw,'// number of cols in matrix SPI\n');
fprintf(fw,'int SD1I_ROWS = %d; ',                  data_CPII.SD1I_ROWS);
fprintf(fw,'// number of rows in matrix SD1I\n');
fprintf(fw,'int SD1I_COLS = %d; ',                  data_CPII.SD1I_COLS);
fprintf(fw,'// number of cols in matrix SD1I\n');
fprintf(fw,'int SD2I_ROWS = %d; ',                  data_CPII.SD2I_ROWS);
fprintf(fw,'// number of rows in matrix SD2I\n');
fprintf(fw,'int SD2I_COLS = %d; ',                  data_CPII.SD2I_COLS);
fprintf(fw,'// number of cols in matrix SD2I\n\n');

% export 1D arrays to a header file
generate_vector(filename, data_CPII.x,              'x',              'double',  data_CPII.x_LEN);
generate_vector(filename, data_CPII.x_prev,         'x_prev',         'double',  data_CPII.x_LEN);
generate_vector(filename, data_CPII.xbar,           'xbar',           'double',  data_CPII.x_LEN);
generate_vector(filename, data_CPII.xbar_prev,      'xbar_prev',      'double',  data_CPII.x_LEN);
generate_vector(filename, data_CPII.nu,             'nu',             'double',  data_CPII.nu_LEN);
generate_vector(filename, data_CPII.nu_prev,        'nu_prev',        'double',  data_CPII.nu_LEN);
generate_vector(filename, data_CPII.p,              'p',              'double',  data_CPII.p_LEN);
generate_vector(filename, data_CPII.p_prev,         'p_prev',         'double',  data_CPII.p_LEN);
generate_vector(filename, data_CPII.p_nu,           'p_nu',           'double',  data_CPII.p_nu_LEN);
generate_vector(filename, data_CPII.p_nu_prev,      'p_nu_prev',      'double',  data_CPII.p_nu_LEN);
generate_vector(filename, data_CPII.b,              'b',              'double',  data_CPII.b_LEN);
generate_vector(filename, data_CPII.f,              'f',              'double',  data_CPII.f_LEN);
generate_vector(filename, data_CPII.c,              'c',              'double',  data_CPII.c_LEN);
generate_vector(filename, data_CPII.l,              'l',              'double',  data_CPII.l_LEN);
generate_vector(filename, data_CPII.diagQ,          'diagQ',          'double',  data_CPII.diagQ_LEN);
generate_vector(filename, data_CPII.new_prox_beg,   'new_prox_beg',   'int',     data_CPII.nProx);
generate_vector(filename, data_CPII.new_prox_end,   'new_prox_end',   'int',     data_CPII.nProx);
generate_vector(filename, data_CPII.len_of_vectors, 'len_of_vectors', 'int',     data_CPII.nProx);

% export matrices to a header file
generate_matrix(filename, data_CPII.A,        'A',        'double',  data_CPII.A_ROWS,        data_CPII.A_COLS);
generate_matrix(filename, data_CPII.K,        'K',        'double',  data_CPII.K_ROWS,        data_CPII.K_COLS);
generate_matrix(filename, data_CPII.KT,       'KT',       'double',  data_CPII.KT_ROWS,       data_CPII.KT_COLS);
generate_matrix(filename, data_CPII.L,        'L',        'double',  data_CPII.L_ROWS,        data_CPII.L_COLS);
generate_matrix(filename, data_CPII.SP,       'SP',       'double',  data_CPII.SP_ROWS,       data_CPII.SP_COLS);
generate_matrix(filename, data_CPII.SD1,      'SD1',      'double',  data_CPII.SD1_ROWS,      data_CPII.SD1_COLS);
generate_matrix(filename, data_CPII.SD2,      'SD2',      'double',  data_CPII.SD2_ROWS,      data_CPII.SD2_COLS);

% export matrice to a header file - indices of nonzero elements
generate_matrix(filename, data_CPII.AI,        'AI',        'int',   data_CPII.AI_ROWS,        data_CPII.AI_COLS);
generate_matrix(filename, data_CPII.KI,        'KI',        'int',   data_CPII.KI_ROWS,        data_CPII.KI_COLS);
generate_matrix(filename, data_CPII.KTI,       'KTI',       'int',   data_CPII.KTI_ROWS,       data_CPII.KTI_COLS);
generate_matrix(filename, data_CPII.LI,        'LI',        'int',   data_CPII.LI_ROWS,        data_CPII.LI_COLS);
generate_matrix(filename, data_CPII.SPI,       'SPI',       'int',   data_CPII.SPI_ROWS,       data_CPII.SPI_COLS);
generate_matrix(filename, data_CPII.SD1I,      'SD1I',      'int',   data_CPII.SD1I_ROWS,      data_CPII.SD1I_COLS);
generate_matrix(filename, data_CPII.SD2I,      'SD2I',      'int',   data_CPII.SD2I_ROWS,      data_CPII.SD2I_COLS);

% generate array of proxes and norms
generate_array_of_prox_types(filename, data_CPII.proxConj_type, data_CPII.nProx);
generate_array_of_norm_types(filename, data_CPII.dualNormType,  data_CPII.nProx);

% generate 1D arrays containing normTypesm, weights and constants to perform projections
generate_vector(filename, data_CPII.weight,    'weights',   'double',  data_CPII.nProx);
generate_vector(filename, data_CPII.constants, 'constants', 'double',  data_CPII.nProx);

% preprocessor directive ends
fw = fopen(filename, 'a');
fprintf(fw, '\n');
fprintf(fw, '#endif');

% close the file handler, otherwise std::exception is thrown
fclose(fw);