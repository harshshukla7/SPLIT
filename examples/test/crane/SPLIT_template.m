make_template('type','PL','project_name','split_project', 'input', ' state0:8:float', 'input', ' primal0:78:float', 'input', ' dual0:28:float', 'input', ' tol_iterates:4:float', 'output', ' primal:78:float', 'output', ' dual:28:float', 'output', ' aux_primal:28:float', 'output', ' iterates:4:float');
 
current_path = pwd; 
 
dest_path_c = strcat(current_path, '/ip_design/src/user_fama.cpp');
dest_path_h = strcat(current_path, '/ip_design/src/user_ama.h');
copyfile('/home/hs/Repository/SPLIT_github/split/coder/src/FPGA/fama_fpga.c', dest_path_c);
copyfile('/home/hs/Repository/SPLIT_github/split/coder/src/FPGA/ama_fpga.h', dest_path_h);
 
 
dest_path_c = strcat(current_path, '/ip_design/src/user_matrix_ops.cpp');
dest_path_h = strcat(current_path, '/ip_design/src/user_matrix_ops.h');
copyfile('/home/hs/Repository/SPLIT_github/split/coder/src/FPGA/matrix_ops_fpga.c', dest_path_c);
copyfile('/home/hs/Repository/SPLIT_github/split/coder/src/FPGA/matrix_ops_fpga.h', dest_path_h);
 
 
source_c = strcat(current_path, '/user_probData.c');
source_h = strcat(current_path, '/user_probData.h');
dest_path_c = strcat(current_path, '/ip_design/src/user_probData.cpp');
dest_path_h = strcat(current_path, '/ip_design/src/user_probData.h');
copyfile(source_c, dest_path_c);
copyfile(source_h, dest_path_h);
 
 
dest_path_c = strcat(current_path, '/ip_design/src/foo_user.cpp');
copyfile('/home/hs/Repository/SPLIT_github/split/coder/src/FPGA/foo_user_fama_orig_fpga.cpp', dest_path_c);
 
 
source_c = strcat(current_path, '/user_mv_mult.cpp');
source_h = strcat(current_path, '/user_mv_mult.h');
dest_path_c = strcat(current_path, '/ip_design/src/user_mv_mult.cpp');
dest_path_h = strcat(current_path, '/ip_design/src/user_mv_mult.h');
copyfile(source_c, dest_path_c);
copyfile(source_h, dest_path_h);
 
 
dest_m = strcat(current_path, '/ip_design/src/test_HIL.m'); 
copyfile('/home/hs/Repository/SPLIT_github/split/FPGA_codegen/test_HIL_orig.m', dest_m); 
 
source_dat = strcat(current_path, '/input_protosplit.mat');
desti_dat = strcat(current_path, '/ip_design/src/input_protosplit.mat');
copyfile(source_dat, desti_dat); 
 
ip_design_build('project_name','split_project','fclk', 100);
 
%%% ip_design_build_debug('project_name','split_project'); 
 
ip_prototype_build('project_name','split_project','board_name','zedboard'); 
 
ip_prototype_load('project_name','split_project','board_name','zedboard','type_eth','udp');
 
ip_prototype_test('project_name','split_project','board_name','zedboard','num_test',1);
 
load('ip_design/src/output_protosplit.mat'); 
