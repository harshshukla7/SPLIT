make_template('type','PL','project_name','split_project', 'input', ' state0:3:float', 'input', ' primal0:5:float', 'input', ' dual0:5:float', 'input', ' tol_iterates:4:float', 'output', ' primal:5:float', 'output', ' dual:5:float', 'output', ' aux_primal:5:float', 'output', ' iterates:4:float');
 
current_path = pwd; 
 
dest_path_c = strcat(current_path, '/ip_design/src/user_admm.cpp');
dest_path_h = strcat(current_path, '/ip_design/src/user_admm.h');
copyfile('/home/hs/Repository/SPLIT_github/split/coder/src/FPGA/admm_fpga.c', dest_path_c);
copyfile('/home/hs/Repository/SPLIT_github/split/coder/src/FPGA/admm_fpga.h', dest_path_h);
 
 
dest_path_c = strcat(current_path, '/ip_design/src/user_matrix_ops.cpp');
dest_path_h = strcat(current_path, '/ip_design/src/user_matrix_ops.h');
copyfile('/home/hs/Repository/SPLIT_github/split/coder/src/FPGA/matrix_ops_fpga.c', dest_path_c);
copyfile('/home/hs/Repository/SPLIT_github/split/coder/src/FPGA/matrix_ops_fpga.h', dest_path_h);
 
 
dest_path_c = strcat(current_path, '/ip_design/src/user_prob_data.cpp');
dest_path_h = strcat(current_path, '/ip_design/src/user_prob_data.h');
copyfile(current_path, dest_path_c);
copyfile(current_path, dest_path_h);
 
 
dest_path_c = strcat(current_path, '/ip_design/src/foo_user.cpp');
copyfile('/home/hs/Repository/SPLIT_github/split/coder/src/FPGA/foo_user_orig.cpp', dest_path_c);
 
 
ip_design_build('project_name','split_project','fclk', 100);
 
%%% ip_design_build_debug('project_name','split_project'); 
 
ip_prototype_build('project_name','split_project','board_name','zedboard'); 
 
ip_prototype_load('project_name','split_project','board_name','zedboard','type_eth','udp');
 
ip_prototype_test('project_name','split_project','board_name','zedboard','num_test',1);
 
