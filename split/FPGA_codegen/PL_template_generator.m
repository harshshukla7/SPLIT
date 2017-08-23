function [output] = PL_template_generator(nParam, nPrimal, nDual, settings)
%% Can have the following arguments
% nParam: number of parameters
% nPrimal: number of primal
% nDual: number of dual
% 1. settings.hardware: 'PL' for pure FPGA and 'soc' for system on chip,
% default: 'PL'
% 2. settings.data_type: 'float' or 'fixed', default: 'float'
% 3. settings.proj_name: any string, default name: 'split_project'
% 4. settings.target_platform: check protoip which platforms are supported,
% default platform: 'zedboard'
% 5. settings.target_frequency: string value of the target frequency (in MHz),
% default:100MHz
% 6. settings.number_test: character value, default:1
% 7. settings.algorithm: character value, default: admm_fpga.c
%    other options: (i)admm (ii)fadmm, (iii) fama, (iv)
%    pda.c (v) fpda


%% hardware

default_hardware = 'PL';
settings.hardware = default_hardware;

%% What kind of data type
default_data_type = 'float';

if ~isfield(settings, 'data_type'), settings.data_type = default_data_type; end
if (strcmp(settings.data_type, 'fixed') == 1)
    
    default_integ_bits = '15';
    default_fract_bits = '15';
    
    if ~isfield(settings, 'integ_bits'), settings.integ_bits = default_integ_bits; end
    if ~isfield(settings, 'fract_bits'), settings.fract_bits = default_fract_bits; end
    
    
end

%% project_name

default_proj_name = 'split_project';
if ~isfield(settings, 'proj_name'), settings.proj_name = default_proj_name; end

%% target platform

default_target_platform = 'zedboard';
if ~isfield(settings, 'target_platform'), settings.target_platform = default_target_platform; end

%% target clock frequency

default_target_frequency = '100';
if ~isfield(settings, 'target_frequency'), settings.target_frequency = default_target_frequency;
else
    if (ischar(settings, 'target_frequency') == 0)
        
        error('please provide target frequency as a character, not a number or in any other format')
    end
end

%% number of tests

default_num_test = '1';
if ~isfield(settings, 'num_tests'), settings.num_tests = default_num_test; end


%% set up for copy paste file

default_algorithm = 'admm';

if ~isfield(settings, 'algorithm'), settings.algorithm = default_algorithm; end

algo = settings.algorithm;

mat_c = which('matrix_ops_fpga.c');
mat_h = which('matrix_ops_fpga.h');
mat_concat_c = '/ip_design/src/user_matrix_ops.cpp';
mat_concat_h = '/ip_design/src/user_matrix_ops.h';

pd_c = '/user_probData.c';
pd_h = '/user_probData.h';
pd_concat_c = '/ip_design/src/user_probData.cpp';
pd_concat_h = '/ip_design/src/user_probData.h';


if ( strcmp(algo, 'admm') == 1)
    
    algo_c = which('admm_fpga.c');
    algo_h = which('admm_fpga.h');
    path_concat_c = '/ip_design/src/user_admm.cpp';
    path_concat_h = '/ip_design/src/user_admm.h';
    foo_concat_c = '/ip_design/src/foo_user.cpp';
    foo_algo_c = which('foo_user_admm_orig_fpga.cpp');
   
    
elseif ( strcmp(algo, 'fadmm') == 1)
    
    algo_c = which('fadmm_fpga.c');
    algo_h = which('admm_fpga.h');
    path_concat_c = '/ip_design/src/user_fadmm.cpp';
    path_concat_h = '/ip_design/src/user_admm.h';
   foo_concat_c = '/ip_design/src/foo_user.cpp';
foo_algo_c = which('foo_user_fadmm_orig_fpga.cpp');
    
elseif ( strcmp(algo, 'fama') == 1)
    
    algo_c = which('fama_fpga.c');
    algo_h = which('ama_fpga.h');
    path_concat_c = '/ip_design/src/user_fama.cpp';
    path_concat_h = '/ip_design/src/user_ama.h';
   foo_concat_c = '/ip_design/src/foo_user.cpp';
foo_algo_c = which('foo_user_fama_orig_fpga.cpp');
    
elseif ( strcmp(algo, 'pda') == 1)
    
    algo_c = which('pda_fpga.c');
    algo_h = which('pda_fpga.h');
    path_concat_c = '/ip_design/src/user_pda.cpp';
    path_concat_h = '/ip_design/src/user_pda.h';
    foo_concat_c = '/ip_design/src/foo_user.cpp';
foo_algo_c = which('foo_user_pda_orig_fpga.cpp');

elseif ( strcmp(algo, 'fpda') == 1)
    
    algo_c = which('fpda_fpga.c');
    algo_h = which('pda_fpga.h');
    path_concat_c = '/ip_design/src/user_fpda.cpp';
    path_concat_h = '/ip_design/src/user_pda.h';
    foo_concat_c = '/ip_design/src/foo_user.cpp';
foo_algo_c = which('foo_user_fpda_orig_fpga.cpp');
    
else
    
    error('unknown algorithm');
    
    
end



%% setting up for code generartion

np = nPrimal;
nd = nDual;
nx = nParam;
hard = settings.hardware;
proj = settings.proj_name;
clk = settings.target_frequency;
board = settings.target_platform;


if strcmp(settings.data_type, 'float')
    
    data_t = 'float';
    
elseif strcmp(settings.data_type, 'fixed')
    data_t = strcat('fix:',settings.integ_bits, ':', settings.fract_bits);
    
else
    error('unknown data type');
end


%% code generation of template file

fileID = fopen('SPLIT_template.m','w');

fprintf(fileID, 'make_template(''type'',''%s'',''project_name'',''%s'', ''input'', '' state0:%d:%s'', ''input'', '' primal0:%d:%s'', ''input'', '' dual0:%d:%s'', ''input'', '' tol_iterates:4:%s'', ''output'', '' primal:%d:%s'', ''output'', '' dual:%d:%s'', ''output'', '' aux_primal:%d:%s'', ''output'', '' iterates:4:%s'');\n \n', hard, proj, nx, data_t, np, data_t, nd, data_t, data_t,  np, data_t, nd, data_t, nd, data_t, data_t  );

%%%%% copy files before start synthesising

fprintf(fileID, 'current_path = pwd; \n \n');

%%%% copy algorithm
fprintf(fileID, 'dest_path_c = strcat(current_path, ''%s'');\n', path_concat_c);
fprintf(fileID, 'dest_path_h = strcat(current_path, ''%s'');\n', path_concat_h);
fprintf(fileID, 'copyfile(''%s'', dest_path_c);\n', algo_c );
fprintf(fileID, 'copyfile(''%s'', dest_path_h);\n \n \n', algo_h );

%%%% copy matrix operations
fprintf(fileID, 'dest_path_c = strcat(current_path, ''%s'');\n', mat_concat_c);
fprintf(fileID, 'dest_path_h = strcat(current_path, ''%s'');\n', mat_concat_h);
fprintf(fileID, 'copyfile(''%s'', dest_path_c);\n', mat_c );
fprintf(fileID, 'copyfile(''%s'', dest_path_h);\n \n \n', mat_h );


%%% problem data copy paste
fprintf(fileID, 'source_c = strcat(current_path, ''%s'');\n', pd_c);
fprintf(fileID, 'source_h = strcat(current_path, ''%s'');\n', pd_h);
fprintf(fileID, 'dest_path_c = strcat(current_path, ''%s'');\n', pd_concat_c);
fprintf(fileID, 'dest_path_h = strcat(current_path, ''%s'');\n', pd_concat_h);
fprintf(fileID, 'copyfile(source_c, dest_path_c);\n' );
fprintf(fileID, 'copyfile(source_h, dest_path_h);\n \n \n' );

%%% foo user copy paste
fprintf(fileID, 'dest_path_c = strcat(current_path, ''%s'');\n', foo_concat_c);
fprintf(fileID, 'copyfile(''%s'', dest_path_c);\n \n \n', foo_algo_c );


%%%%% copy ends

fprintf(fileID, 'ip_design_build(''project_name'',''%s'',''fclk'', %s);\n \n', proj, clk);

fprintf(fileID, '%%%%%% ip_design_build_debug(''project_name'',''%s''); \n \n', proj);

fprintf(fileID, 'ip_prototype_build(''project_name'',''%s'',''board_name'',''%s''); \n \n',proj, board);

fprintf(fileID, 'ip_prototype_load(''project_name'',''%s'',''board_name'',''%s'',''type_eth'',''udp'');\n \n', proj, board);
fprintf(fileID, 'ip_prototype_test(''project_name'',''%s'',''board_name'',''%s'',''num_test'',%s);\n \n',proj, board, settings.num_tests);
fclose(fileID);
output = 1;
%% delete all the folllowings
% if (isfield (settings))
% fprintf(fileID,'#define Size_Rows %d\n', Size_rows);
% fprintf(fileID,'#define Size_Columns %d\n', Size_columns);
% fprintf(fileID,'#define Num_PARAL %d\n', Num_Par);
% fprintf(fileID,'#define PART_Size %d\n', Part_size);
% fprintf(fileID,'#define REM_PART_Size %d\n', Rem_Part_size);
% fprintf(fileID,'#define ACC_Size %d\n\n', number_y_local);
%
% fprintf(fileID,strcat('void mv_mult(',data_t,' y_out[SIZE],',data_t,' x_in[SIZE]);\n'));
% fclose(fileID);


end