%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% icl::protoip
% Author: asuardi <https://github.com/asuardi>
% Date: November - 2014
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [pl_vec_out_out_int] = foo_user(project_name,pl_vec_in_in_int, pl_mat_in_in_int)


	% load project configuration parameters: input and output vectors (name, size, type, NUM_TEST, TYPE_TEST)
	load_configuration_parameters(project_name);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% compute with Matlab and save in a file simulation results pl_vec_out_out_int
	for i=1:PL_VEC_OUT_OUT_LENGTH
		pl_vec_out_out_int(i)=0;
		for i_pl_vec_in = 1:PL_VEC_IN_IN_LENGTH
			pl_vec_out_out_int(i)=pl_vec_out_out_int(i) + pl_vec_in_in_int(i_pl_vec_in);
		end
		for i_pl_mat_in = 1:PL_MAT_IN_IN_LENGTH
			pl_vec_out_out_int(i)=pl_vec_out_out_int(i) + pl_mat_in_in_int(i_pl_mat_in);
		end
	end

end

