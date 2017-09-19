make_template('type','SOC','project_name','split_project', 'input', ' pl_vec_in:5:float', 'input', ' pl_mat_in:25:float', 'output', ' pl_vec_out:5:float', 'soc_input', ' state0:5:float', 'soc_input', ' primal0:5:float', 'soc_input', ' dual0:5:float', 'soc_input', ' tol_iterates:4:float', 'soc_output', ' primal:5:float', 'soc_output', ' dual:5:float', 'soc_output', ' aux_primal:5:float', 'soc_output', ' iterates:4:float');
 
current_path = pwd; 
 
dest_path_c = strcat(current_path, '/soc_prototype/src/user_admm.c');
dest_path_h = strcat(current_path, '/soc_prototype/src/user_admm.h');
copyfile('/home/hs/Repository/SPLIT_github/split/coder/src/SoC/admm_soc.c', dest_path_c);
copyfile('/home/hs/Repository/SPLIT_github/split/coder/src/SoC/admm_soc.h', dest_path_h);
 
 
dest_path_c = strcat(current_path, '/soc_prototype/src/user_matrix_ops.c');
dest_path_h = strcat(current_path, '/soc_prototype/src/user_matrix_ops.h');
copyfile('/home/hs/Repository/SPLIT_github/split/coder/src/FPGA/matrix_ops_fpga.c', dest_path_c);
copyfile('/home/hs/Repository/SPLIT_github/split/coder/src/FPGA/matrix_ops_fpga.h', dest_path_h);
 
 
dest_path_c = strcat(current_path, '/soc_prototype/src/user_prob_data.c');
dest_path_h = strcat(current_path, '/soc_prototype/src/user_prob_data.h');
copyfile(current_path, dest_path_c);
copyfile(current_path, dest_path_h);
 
 
dest_path_c = strcat(current_path, '/ip_design/src/user_mat_vec.cpp');
dest_path_h = strcat(current_path, '/ip_design/src/user_mat_vec.h');
copyfile(current_path, dest_path_c);
copyfile(current_path, dest_path_h);
 
 
dest_path_c = strcat(current_path, '/ip_design/src/foo_user.cpp');
copyfile('/home/hs/Repository/SPLIT_github/split/coder/src/SoC/foo_user_orig_soc.cpp', dest_path_c);
 
 
dest_path_c = strcat(current_path, '/soc_prototype/src/soc_user.c');
copyfile('/home/hs/Repository/SPLIT_github/split/coder/src/SoC/soc_user_admm_orig.c', dest_path_c);
 
 
ip_design_build('project_name','split_project','fclk', 100);
 
%%% ip_design_build_debug('project_name','split_project'); 
 
ip_prototype_build('project_name','split_project','board_name','zedboard'); 
 
soc_prototype_load('project_name','split_project','board_name','zedboard','type_eth','udp');
 
%%% soc_prototype_load_debug('project_name','split_project','board_name','zedboard');
 
soc_prototype_test('project_name','split_project','board_name','zedboard','num_test',1);
 
