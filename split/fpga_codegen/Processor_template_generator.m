function [output] = Processor_template_generator(nParam, nPrimal, nDual, settings)
%% Can have the following arguments
% nstates: number of states
% nPrimal: number of primal
% nDual: number of dual
% 1. settings.hardware: 'SOC' 
% 2. settings.data_type: 'float' or 'fixed', default: 'float'
% 3. settings.proj_name: any string, default name: 'split_project'
% 4. settings.target_platform: check protoip which platforms are supported,
% default platform: 'zedboard'
% 5. settings.target_frequency: string value of the target frequency (in MHz),
% default:100MHz
% 6. settings.number_test: character value, default:1
% 7. settings.x0
% 8. settings.p0
% 9. settings.d0
% 10.settings.resi_tol  
% 11.settings.algorithm: character value, default: admm_fpga.c
%    other options: (i)admm (ii)fadmm, (iii) fama, (iv)
%    pda.c (v) fpda




%% hardware

default_hardware = 'SOC'; %Embedded Processor
settings.hardware = default_hardware;


%% What kind of data type
default_data_type = 'float';

if ~isfield(settings, 'data_type'), settings.data_type = default_data_type; end
if (strcmp(settings.data_type, 'fixed') == 1)
    
    error('Fixed point is currently not supported for Embedded Processors')
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


%% check which inputs variables are passed else assign some values!

if ~isfield(settings, 'x0'), settings.x0 = zeros(1,nParam); end
if ~isfield(settings, 'p0'), settings.p0 = zeros(1,nPrimal); end
if ~isfield(settings, 'd0'), settings.d0 = zeros(1,nDual); end
if ~isfield(settings, 'resi_tol'), settings.resi_tol = [1e-2, 1e-2, 200, 20]; end

x0 = settings.x0;
p0 = settings.p0;
d0 = settings.d0;
tol_itr0 = settings.resi_tol;

save('input_protosplit.mat', 'x0', 'p0', 'd0', 'tol_itr0');


%% set up for copy paste file


%%% matrix operations
mat_c = which('matrix_ops_ep.c');
mat_h = which('matrix_ops_ep.h');
mat_concat_c = '/soc_prototype/src/user_matrix_ops.c';
mat_concat_h = '/soc_prototype/src/user_matrix_ops.h';

%%% problem data
pd_c = '/user_probData.c';
pd_h = '/user_probData.h';
pd_concat_c = '/soc_prototype/src/user_probData.c';
pd_concat_h = '/soc_prototype/src/user_probData.h';


%%% foo user for ip design
foo_concat_c = '/ip_design/src/foo_user.cpp';
foo_algo_c = which('foo_user_orig_ep.cpp');


%%% foo user for soc prototype
foo_soc_concat_c = '/soc_prototype/src/user_foo_data.h';
foo_soc_algo_c = which('user_foo_data_ep_protoype.h');

%%% data declaration for soc
data_decl_source_soc = which('user_foo_data_ep_protoype.h');
data_decl_concat_soc = ('/soc_prototype/src/user_foo_data.h');


%%% test_HIL for soc
th_source_soc = which('test_HIL_orig_ep.m');
th_concat_soc = ('/soc_prototype/src/test_HIL.m');

%%% test_HIL for fpga
th_source_fpga = which('test_HIL_orig_ep_fpga.m');
th_concat_fpga = ('/ip_design/src/test_HIL.m');

%%% ldl.c for suitesparse
ldlc_source_soc = which('user_ldl_orig.c');
ldlc_concat_soc = ('/soc_prototype/src/user_ldl.c');


%%% ldl.h for suitesparse
ldlh_source_soc = which('user_ldl_orig.h');
ldlh_concat_soc = ('/soc_prototype/src/user_ldl.h');


%%% suitesparse config
ss_config_source_soc = which('SuiteSparse_config_orig.h');
ss_config_concat_soc = ('/soc_prototype/src/user_suiteSparse_config.h');


%%% copy paste input data file
input_data_source_concat = '/input_protosplit.mat';
input_data_dest_concat = '/soc_prototype/src/input_protosplit.mat';


default_algorithm = 'admm';

if ~isfield(settings, 'algorithm'), settings.algorithm = default_algorithm; end

algo = settings.algorithm;


if ( strcmp(algo, 'admm') == 1)
    
    algo_c = which('admm_ep.c');
    algo_h = which('admm_ep.h');
    path_concat_c = '/soc_prototype/src/user_admm.c';
    path_concat_h = '/soc_prototype/src/user_admm.h';

    soc_concat_c = '/soc_prototype/src/soc_user.c';
    soc_algo_c = which('soc_user_admm_orig_ep.c');

    
    
elseif ( strcmp(algo, 'fadmm') == 1)
    
    algo_c = which('fadmm_ep.c');
    algo_h = which('admm_ep.h');
    path_concat_c = '/soc_prototype/src/user_fadmm.c';
    path_concat_h = '/soc_prototype/src/user_admm.h';

    soc_concat_c = '/soc_prototype/src/soc_user.c';
    soc_algo_c = which('soc_user_fadmm_orig_ep.c');

    
elseif ( strcmp(algo, 'fama') == 1)
    
    algo_c = which('fama_ep.c');
    algo_h = which('ama_ep.h');
    path_concat_c = '/soc_prototype/src/user_fama.c';
    path_concat_h = '/soc_prototype/src/user_ama.h';

    soc_concat_c = '/soc_prototype/src/soc_user.c';
    soc_algo_c = which('soc_user_fama_orig_ep.c');

    
elseif ( strcmp(algo, 'pda') == 1)
    
    algo_c = which('pda_ep.c');
    algo_h = which('pda_ep.h');
    path_concat_c = '/soc_prototype/src/user_pda.c';
    path_concat_h = '/soc_prototype/src/user_pda.h';

    soc_concat_c = '/soc_prototype/src/soc_user.c';
    soc_algo_c = which('soc_user_pda_orig_ep.c');

    
elseif ( strcmp(algo, 'fpda') == 1)
    
    algo_c = which('fpda_ep.c');
    algo_h = which('pda_ep.h');
    path_concat_c = '/soc_prototype/src/user_fpda.c';
    path_concat_h = '/soc_prototype/src/user_pda.h';

    soc_concat_c = '/soc_prototype/src/soc_user.c';
    soc_algo_c = which('soc_user_fpda_orig_ep.c');

    
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

fprintf(fileID, 'make_template(''type'',''%s'',''project_name'',''%s'', ''input'', '' pl_vec_in:1:float'', ''output'', '' pl_vec_out:1:float'', ''soc_input'', '' state0:%d:%s'', ''soc_input'', '' primal0:%d:%s'', ''soc_input'', '' dual0:%d:%s'', ''soc_input'', '' tol_iterates:4:%s'', ''soc_output'', '' primal:%d:%s'', ''soc_output'', '' dual:%d:%s'', ''soc_output'', '' aux_primal:%d:%s'', ''soc_output'', '' iterates:4:%s'');\n \n', hard, proj,  nx, data_t, np, data_t, nd, data_t, data_t,  np, data_t, nd, data_t, nd, data_t, data_t  );

%%% following line was for sending matrix as well
%fprintf(fileID, 'make_template(''type'',''%s'',''project_name'',''%s'', ''input'', '' pl_vec_in:%d:%s'', ''input'', '' pl_mat_in:%d:%s'', ''output'', '' pl_vec_out:%d:%s'', ''soc_input'', '' state0:%d:%s'', ''soc_input'', '' primal0:%d:%s'', ''soc_input'', '' dual0:%d:%s'', ''soc_input'', '' tol_iterates:4:%s'', ''soc_output'', '' primal:%d:%s'', ''soc_output'', '' dual:%d:%s'', ''soc_output'', '' aux_primal:%d:%s'', ''soc_output'', '' iterates:4:%s'');\n \n', hard, proj, input_dim, data_t, (mat_row*mat_col), data_t, output_dim, data_t, nx, data_t, np, data_t, nd, data_t, data_t,  np, data_t, nd, data_t, nd, data_t, data_t  );

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


%%%data delcaration copy paste
fprintf(fileID, 'dest_path_c = strcat(current_path, ''%s'');\n', data_decl_concat_soc);
fprintf(fileID, 'copyfile(''%s'', dest_path_c);\n \n \n', data_decl_source_soc );


%%% foo user for ip design copy paste
fprintf(fileID, 'dest_path_c = strcat(current_path, ''%s'');\n', foo_concat_c);
fprintf(fileID, 'copyfile(''%s'', dest_path_c);\n \n \n', foo_algo_c );

%%% foo user for soc prototype copy paste
fprintf(fileID, 'dest_path_c = strcat(current_path, ''%s'');\n', foo_soc_concat_c);
fprintf(fileID, 'copyfile(''%s'', dest_path_c);\n \n \n', foo_soc_algo_c );


%%% soc user copy paste
fprintf(fileID, 'dest_path_c = strcat(current_path, ''%s'');\n', soc_concat_c);
fprintf(fileID, 'copyfile(''%s'', dest_path_c);\n \n \n', soc_algo_c );

%%% copy paste test_hil for soc_prototype
fprintf(fileID, 'dest_m = strcat(current_path, ''%s''); \n', th_concat_soc);
fprintf(fileID, 'copyfile(''%s'', dest_m); \n \n', th_source_soc);

%%% copy paste test_hil for ip_design
fprintf(fileID, 'dest_m = strcat(current_path, ''%s''); \n', th_concat_fpga);
fprintf(fileID, 'copyfile(''%s'', dest_m); \n \n', th_source_fpga);

%%% copy paste ldl.c
fprintf(fileID, 'dest_m = strcat(current_path, ''%s''); \n', ldlc_concat_soc);
fprintf(fileID, 'copyfile(''%s'', dest_m); \n \n', ldlc_source_soc);

%%% copy paste ldl.h
fprintf(fileID, 'dest_m = strcat(current_path, ''%s''); \n', ldlh_concat_soc);
fprintf(fileID, 'copyfile(''%s'', dest_m); \n \n', ldlh_source_soc);

%%% copy paste suite sparse
fprintf(fileID, 'dest_m = strcat(current_path, ''%s''); \n', ss_config_concat_soc);
fprintf(fileID, 'copyfile(''%s'', dest_m); \n \n', ss_config_source_soc);

%%% copy input data file

fprintf(fileID, 'source_dat = strcat(current_path, ''%s'');\n', input_data_source_concat);
fprintf(fileID, 'desti_dat = strcat(current_path, ''%s'');\n', input_data_dest_concat);

fprintf(fileID, 'copyfile(source_dat, desti_dat); \n \n' );


%%%%% copy ends

fprintf(fileID, 'ip_design_build(''project_name'',''%s'',''fclk'', %s);\n \n', proj, clk);

fprintf(fileID, '%%%%%% ip_design_build_debug(''project_name'',''%s''); \n \n', proj);

fprintf(fileID, 'ip_prototype_build(''project_name'',''%s'',''board_name'',''%s''); \n \n',proj, board);

fprintf(fileID, 'soc_prototype_load(''project_name'',''%s'',''board_name'',''%s'',''type_eth'',''udp'');\n \n', proj, board);
fprintf(fileID, '%%%%%% soc_prototype_load_debug(''project_name'',''%s'',''board_name'',''%s'');\n \n', proj, board);

fprintf(fileID, 'soc_prototype_test(''project_name'',''%s'',''board_name'',''%s'',''num_test'',%s);\n \n',proj, board, settings.num_tests);
fclose(fileID);
output = 1;


while ~exist([pwd filesep 'SPLIT_template.m'], 'file') ; end
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