function export_matrices_parameters_FGM(data_FGM)

% create a file handler and the corresponding header file
filename = 'FGM_params.h';
fw = fopen(filename, 'w');

%% FGM parameters
fprintf(fw,'#define EPS %d ',data_FGM.epsilon);
fprintf(fw,'// relative tolerance\n');
fprintf(fw,'#define STEP %d ',data_FGM.step);
fprintf(fw,'// number of iterations\n');
fprintf(fw,'#define MU %f ',data_FGM.mu);
fprintf(fw,'// convexity parameter\n');
fprintf(fw,'#define L %d ',data_FGM.L);
fprintf(fw,'// Lipschitz constant \n');
fprintf(fw,'#define ALPHA_I %d ',data_FGM.alpha0);
fprintf(fw,'// auxiliary parameter: Nesterov acceleration\n');
fprintf(fw,'#define ALPHA_I_1 %d ',0);
fprintf(fw,'// auxiliary parameter: Nesterov acceleration\n');

%% size of matrices, occuring in the fast gradient method

% number of rows and cols in matrix H_g
fprintf(fw,'#define HG_ROWS %d ', size(data_FGM.H_g, 1));
fprintf(fw,'// number of rows of matrix H_g\n');
fprintf(fw,'#define HG_COLS %d ', size(data_FGM.H_g, 2));
fprintf(fw,'// number of cols of matrix H_g\n');

% number of rows and cols in matrix H_g
fprintf(fw,'#define M1G_ROWS %d ', size(data_FGM.M1_g, 2));
fprintf(fw,'// number of rows of matrix M1_g\n');
fprintf(fw,'#define M1G_COLS %d ', size(data_FGM.M1_g, 1));
fprintf(fw,'// number of cols of matrix M1_g\n');

% number of rows and cols in matrix H_g
fprintf(fw,'#define M2G_ROWS %d ', size(data_FGM.M2_g, 1));
fprintf(fw,'// number of rows of matrix M2_g\n');
fprintf(fw,'#define M2G_COLS %d ', size(data_FGM.M2_g, 2));
fprintf(fw,'// number of cols of matrix M2_g\n\n');


%% export matrices

% matrix H_g
generate_matrix(filename, data_FGM.H_g, 'H_g', size(data_FGM.H_g,1), size(data_FGM.H_g, 2));
% matrix M1_g
generate_matrix(filename, data_FGM.M1_g', 'M1_g', size(data_FGM.M1_g', 1), size(data_FGM.M1_g', 2));
% matrix M2_g
generate_matrix(filename, data_FGM.M2_g, 'M2_g', size(data_FGM.M2_g, 1), size(data_FGM.M2_g, 2));

%% export vectors

% u_g - solution
generate_vector(filename, zeros(data_FGM.horizon * data_FGM.num_inputs, 1), 'u_fg', data_FGM.horizon * data_FGM.num_inputs);
% u_fg - solution from previous step
generate_vector(filename, zeros(data_FGM.horizon * data_FGM.num_inputs, 1), 'u_fg_prev', data_FGM.horizon * data_FGM.num_inputs);
% u_fg_hat - update of beta
generate_vector(filename, zeros(data_FGM.horizon * data_FGM.num_inputs, 1), 'u_fg_hat', data_FGM.horizon * data_FGM.num_inputs);
% u_fg_hat_prev - previus value of u_fg
generate_vector(filename, zeros(data_FGM.horizon * data_FGM.num_inputs, 1), 'u_fg_hat_prev', data_FGM.horizon * data_FGM.num_inputs);
% x0 - initial condition
generate_vector(filename, data_FGM.init_conditions, 'x0', size(data_FGM.init_conditions, 1));

%% box constraints
generate_vector(filename, -ones(data_FGM.horizon * data_FGM.num_inputs, 1), 'umin', data_FGM.horizon * data_FGM.num_inputs);
generate_vector(filename,  ones(data_FGM.horizon * data_FGM.num_inputs, 1), 'umax', data_FGM.horizon * data_FGM.num_inputs);

% ----------------------------------------------------------------------- %
function generate_matrix(file, M, stringM, n, m)

fw = fopen(file,'a');

fprintf(fw,'float %s[%d][%d] = {',stringM, n, m);
for i = 1:n
    fprintf(fw,' {');
    for j = 1:m
        if(j == m)
            fprintf(fw,'%d ', M(i,j));
        else
            fprintf(fw,'%d, ',M(i,j)); 
        end
    end
    if(i == n)
        fprintf(fw,'} ');
    else
        fprintf(fw,'}, ');
    end
end
fprintf(fw,'}; \n');

% ----------------------------------------------------------------------- %
function generate_vector(file, M, stringM, n)

fw = fopen(file,'a');

fprintf(fw,'float %s[%d] = {',stringM, n);
for i = 1:n
    if(i == n)
       fprintf(fw,'%d ' , M(i));
    else
       fprintf(fw,'%d, ', M(i)); 
    end
end
fprintf(fw,'}; \n');