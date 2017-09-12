%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% icl::protoip
% Author: asuardi <https://github.com/asuardi>
% Date: November - 2014
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function test_HIL(project_name)


addpath('../../.metadata');
mex FPGAclientMATLAB.c
load_configuration_parameters(project_name)


% rng('shuffle');

for i=1:NUM_TEST
	tmp_disp_str=strcat('Test number ',num2str(i));
	disp(tmp_disp_str)


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% generate random stimulus vector pl_vec_in_in. (-5<=pl_vec_in_in <=5)
% 	pl_vec_in_in=rand(1,PL_VEC_IN_IN_LENGTH)*10-5;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%save pl_vec_in_in_log
	if (TYPE_TEST==0)
		filename = strcat('../../ip_prototype/test/results/', project_name ,'/pl_vec_in_in_log.dat');
	else
		filename = strcat('../test/results/', project_name ,'/pl_vec_in_in_log.dat');
	end
	fid = fopen(filename, 'a+');
   
	for j=1:length(pl_vec_in_in)
		fprintf(fid, '%2.18f,',pl_vec_in_in(j));
	end
	fprintf(fid, '\n');

	fclose(fid);


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% generate random stimulus vector pl_mat_in_in. (-5<=pl_mat_in_in <=5)
% 	pl_mat_in_in=rand(1,PL_MAT_IN_IN_LENGTH)*10-5;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%save pl_mat_in_in_log
	if (TYPE_TEST==0)
		filename = strcat('../../ip_prototype/test/results/', project_name ,'/pl_mat_in_in_log.dat');
	else
		filename = strcat('../test/results/', project_name ,'/pl_mat_in_in_log.dat');
	end
	fid = fopen(filename, 'a+');
   
	for j=1:length(pl_mat_in_in)
		fprintf(fid, '%2.18f,',pl_mat_in_in(j));
	end
	fprintf(fid, '\n');

	fclose(fid);


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Start Matlab timer
	tic

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% send the stimulus to the FPGA simulation model when IP design test or to FPGA evaluation borad when IP prototype, execute the algorithm and read back the results
	% reset IP
	Packet_type=1; % 1 for reset, 2 for start, 3 for write to IP vector packet_internal_ID, 4 for read from IP vector packet_internal_ID of size packet_output_size
	packet_internal_ID=0;
	packet_output_size=1;
	data_to_send=1;
	FPGAclientMATLAB(data_to_send,Packet_type,packet_internal_ID,packet_output_size);


	% send data to FPGA
	% send pl_vec_in_in
	Packet_type=3; % 1 for reset, 2 for start, 3 for write to IP vector packet_internal_ID, 4 for read from IP vector packet_internal_ID of size packet_output_size
	packet_internal_ID=0;
	packet_output_size=1;
	data_to_send=pl_vec_in_in;
	FPGAclientMATLAB(data_to_send,Packet_type,packet_internal_ID,packet_output_size);

	% send pl_mat_in_in
	Packet_type=3; % 1 for reset, 2 for start, 3 for write to IP vector packet_internal_ID, 4 for read from IP vector packet_internal_ID of size packet_output_size
	packet_internal_ID=1;
	packet_output_size=1;
	data_to_send=pl_mat_in_in;
	FPGAclientMATLAB(data_to_send,Packet_type,packet_internal_ID,packet_output_size);


	% start FPGA
	Packet_type=2; % 1 for reset, 2 for start, 3 for write to IP vector packet_internal_ID, 4 for read from IP vector packet_internal_ID of size packet_output_size
	packet_internal_ID=0;
	packet_output_size=1;
	data_to_send=0;
	FPGAclientMATLAB(data_to_send,Packet_type,packet_internal_ID,packet_output_size);


	% read data from FPGA
	% read fpga_pl_vec_out_out
	Packet_type=4; % 1 for reset, 2 for start, 3 for write to IP vector packet_internal_ID, 4 for read from IP vector packet_internal_ID of size packet_output_size
	packet_internal_ID=0;
	packet_output_size=PL_VEC_OUT_OUT_LENGTH;
	data_to_send=0;
	[output_FPGA, time_IP] = FPGAclientMATLAB(data_to_send,Packet_type,packet_internal_ID,packet_output_size);
	fpga_pl_vec_out_out=output_FPGA;
	% Stop Matlab timer
	time_matlab=toc;
	time_communication=time_matlab-time_IP;


	%save fpga_pl_vec_out_out_log.dat
	if (TYPE_TEST==0)
		filename = strcat('../../ip_prototype/test/results/', project_name ,'/fpga_pl_vec_out_out_log.dat');
	else
		filename = strcat('../test/results/', project_name ,'/fpga_pl_vec_out_out_log.dat');
	end
	fid = fopen(filename, 'a+');
   
	for j=1:length(fpga_pl_vec_out_out)
		fprintf(fid, '%2.18f,',fpga_pl_vec_out_out(j));
	end
	fprintf(fid, '\n');

	fclose(fid);



	%save fpga_time_log.dat
	if (TYPE_TEST==0)
		filename = strcat('../../ip_prototype/test/results/', project_name ,'/fpga_time_log.dat');
	else
		filename = strcat('../test/results/', project_name ,'/fpga_time_log.dat');
	end
	fid = fopen(filename, 'a+');
   
	fprintf(fid, '%2.18f, %2.18f \n',time_IP, time_communication);

	fclose(fid);


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% compute with Matlab and save in a file simulation results
% 	[matlab_pl_vec_out_out] = foo_user(project_name,pl_vec_in_in, pl_mat_in_in);


	%save matlab_pl_vec_out_out_log
	if (TYPE_TEST==0)
		filename = strcat('../../ip_prototype/test/results/', project_name ,'/matlab_pl_vec_out_out_log.dat');
	else
		filename = strcat('../test/results/', project_name ,'/matlab_pl_vec_out_out_log.dat');
	end
	fid = fopen(filename, 'a+');
   
	for j=1:length(matlab_pl_vec_out_out)
		fprintf(fid, '%2.18f,',matlab_pl_vec_out_out(j));
	end
	fprintf(fid, '\n');

	fclose(fid);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%write a dummy file to tell tcl script to continue with the execution

filename = strcat('_locked');
fid = fopen(filename, 'w');
fprintf(fid, 'locked write\n');
fclose(fid);

if strcmp(TYPE_DESIGN_FLOW,'vivado')
	quit;
end

end
