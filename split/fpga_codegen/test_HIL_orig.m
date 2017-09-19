%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% icl::protoip
% Author: asuardi <https://github.com/asuardi>
% Date: November - 2014
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function test_HIL(project_name)


addpath('../../.metadata');
mex FPGAclientMATLAB.c
load_configuration_parameters(project_name)


%%%%% load data from the workspace
load('input_protosplit.mat');


% rng('shuffle');

for i=1:NUM_TEST
	tmp_disp_str=strcat('Test number ',num2str(i));
	disp(tmp_disp_str)


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% generate random stimulus vector state0_in. (-5<=state0_in <=5)
	state0_in=x0;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%save state0_in_log
	if (TYPE_TEST==0)
		filename = strcat('../../ip_prototype/test/results/', project_name ,'/state0_in_log.dat');
	else
		filename = strcat('../test/results/', project_name ,'/state0_in_log.dat');
	end
	fid = fopen(filename, 'a+');
   
	for j=1:length(state0_in)
		fprintf(fid, '%2.18f,',state0_in(j));
	end
	fprintf(fid, '\n');

	fclose(fid);


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% generate random stimulus vector primal0_in. (-5<=primal0_in <=5)
	primal0_in=p0;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%save primal0_in_log
	if (TYPE_TEST==0)
		filename = strcat('../../ip_prototype/test/results/', project_name ,'/primal0_in_log.dat');
	else
		filename = strcat('../test/results/', project_name ,'/primal0_in_log.dat');
	end
	fid = fopen(filename, 'a+');
   
	for j=1:length(primal0_in)
		fprintf(fid, '%2.18f,',primal0_in(j));
	end
	fprintf(fid, '\n');

	fclose(fid);


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% generate random stimulus vector dual0_in. (-5<=dual0_in <=5)
	dual0_in=d0;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%save dual0_in_log
	if (TYPE_TEST==0)
		filename = strcat('../../ip_prototype/test/results/', project_name ,'/dual0_in_log.dat');
	else
		filename = strcat('../test/results/', project_name ,'/dual0_in_log.dat');
	end
	fid = fopen(filename, 'a+');
   
	for j=1:length(dual0_in)
		fprintf(fid, '%2.18f,',dual0_in(j));
	end
	fprintf(fid, '\n');

	fclose(fid);


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% generate random stimulus vector tol_iterates_in. (-5<=tol_iterates_in <=5)
	tol_iterates_in=tol_itr0;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%save tol_iterates_in_log
	if (TYPE_TEST==0)
		filename = strcat('../../ip_prototype/test/results/', project_name ,'/tol_iterates_in_log.dat');
	else
		filename = strcat('../test/results/', project_name ,'/tol_iterates_in_log.dat');
	end
	fid = fopen(filename, 'a+');
   
	for j=1:length(tol_iterates_in)
		fprintf(fid, '%2.18f,',tol_iterates_in(j));
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
	% send state0_in
	Packet_type=3; % 1 for reset, 2 for start, 3 for write to IP vector packet_internal_ID, 4 for read from IP vector packet_internal_ID of size packet_output_size
	packet_internal_ID=0;
	packet_output_size=1;
	data_to_send=state0_in;
	FPGAclientMATLAB(data_to_send,Packet_type,packet_internal_ID,packet_output_size);

	% send primal0_in
	Packet_type=3; % 1 for reset, 2 for start, 3 for write to IP vector packet_internal_ID, 4 for read from IP vector packet_internal_ID of size packet_output_size
	packet_internal_ID=1;
	packet_output_size=1;
	data_to_send=primal0_in;
	FPGAclientMATLAB(data_to_send,Packet_type,packet_internal_ID,packet_output_size);

	% send dual0_in
	Packet_type=3; % 1 for reset, 2 for start, 3 for write to IP vector packet_internal_ID, 4 for read from IP vector packet_internal_ID of size packet_output_size
	packet_internal_ID=2;
	packet_output_size=1;
	data_to_send=dual0_in;
	FPGAclientMATLAB(data_to_send,Packet_type,packet_internal_ID,packet_output_size);

	% send tol_iterates_in
	Packet_type=3; % 1 for reset, 2 for start, 3 for write to IP vector packet_internal_ID, 4 for read from IP vector packet_internal_ID of size packet_output_size
	packet_internal_ID=3;
	packet_output_size=1;
	data_to_send=tol_iterates_in;
	FPGAclientMATLAB(data_to_send,Packet_type,packet_internal_ID,packet_output_size);


	% start FPGA
	Packet_type=2; % 1 for reset, 2 for start, 3 for write to IP vector packet_internal_ID, 4 for read from IP vector packet_internal_ID of size packet_output_size
	packet_internal_ID=0;
	packet_output_size=1;
	data_to_send=0;
	FPGAclientMATLAB(data_to_send,Packet_type,packet_internal_ID,packet_output_size);


	% read data from FPGA
	% read fpga_primal_out
	Packet_type=4; % 1 for reset, 2 for start, 3 for write to IP vector packet_internal_ID, 4 for read from IP vector packet_internal_ID of size packet_output_size
	packet_internal_ID=0;
	packet_output_size=PRIMAL_OUT_LENGTH;
	data_to_send=0;
	[output_FPGA, time_IP] = FPGAclientMATLAB(data_to_send,Packet_type,packet_internal_ID,packet_output_size);
	fpga_primal_out=output_FPGA;
	% read fpga_dual_out
	Packet_type=4; % 1 for reset, 2 for start, 3 for write to IP vector packet_internal_ID, 4 for read from IP vector packet_internal_ID of size packet_output_size
	packet_internal_ID=1;
	packet_output_size=DUAL_OUT_LENGTH;
	data_to_send=0;
	[output_FPGA, time_IP] = FPGAclientMATLAB(data_to_send,Packet_type,packet_internal_ID,packet_output_size);
	fpga_dual_out=output_FPGA;
	% read fpga_aux_primal_out
	Packet_type=4; % 1 for reset, 2 for start, 3 for write to IP vector packet_internal_ID, 4 for read from IP vector packet_internal_ID of size packet_output_size
	packet_internal_ID=2;
	packet_output_size=AUX_PRIMAL_OUT_LENGTH;
	data_to_send=0;
	[output_FPGA, time_IP] = FPGAclientMATLAB(data_to_send,Packet_type,packet_internal_ID,packet_output_size);
	fpga_aux_primal_out=output_FPGA;
	% read fpga_iterates_out
	Packet_type=4; % 1 for reset, 2 for start, 3 for write to IP vector packet_internal_ID, 4 for read from IP vector packet_internal_ID of size packet_output_size
	packet_internal_ID=3;
	packet_output_size=ITERATES_OUT_LENGTH;
	data_to_send=0;
	[output_FPGA, time_IP] = FPGAclientMATLAB(data_to_send,Packet_type,packet_internal_ID,packet_output_size);
	fpga_iterates_out=output_FPGA;
	% Stop Matlab timer
	time_matlab=toc;
	time_communication=time_matlab-time_IP;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% save the variable
    save('output_protosplit.mat','fpga_primal_out', 'fpga_dual_out', 'fpga_aux_primal_out', 'fpga_iterates_out')


	%save fpga_primal_out_log.dat
	if (TYPE_TEST==0)
		filename = strcat('../../ip_prototype/test/results/', project_name ,'/fpga_primal_out_log.dat');
	else
		filename = strcat('../test/results/', project_name ,'/fpga_primal_out_log.dat');
	end
	fid = fopen(filename, 'a+');
   
	for j=1:length(fpga_primal_out)
		fprintf(fid, '%2.18f,',fpga_primal_out(j));
	end
	fprintf(fid, '\n');

	fclose(fid);



	%save fpga_dual_out_log.dat
	if (TYPE_TEST==0)
		filename = strcat('../../ip_prototype/test/results/', project_name ,'/fpga_dual_out_log.dat');
	else
		filename = strcat('../test/results/', project_name ,'/fpga_dual_out_log.dat');
	end
	fid = fopen(filename, 'a+');
   
	for j=1:length(fpga_dual_out)
		fprintf(fid, '%2.18f,',fpga_dual_out(j));
	end
	fprintf(fid, '\n');

	fclose(fid);



	%save fpga_aux_primal_out_log.dat
	if (TYPE_TEST==0)
		filename = strcat('../../ip_prototype/test/results/', project_name ,'/fpga_aux_primal_out_log.dat');
	else
		filename = strcat('../test/results/', project_name ,'/fpga_aux_primal_out_log.dat');
	end
	fid = fopen(filename, 'a+');
   
	for j=1:length(fpga_aux_primal_out)
		fprintf(fid, '%2.18f,',fpga_aux_primal_out(j));
	end
	fprintf(fid, '\n');

	fclose(fid);



	%save fpga_iterates_out_log.dat
	if (TYPE_TEST==0)
		filename = strcat('../../ip_prototype/test/results/', project_name ,'/fpga_iterates_out_log.dat');
	else
		filename = strcat('../test/results/', project_name ,'/fpga_iterates_out_log.dat');
	end
	fid = fopen(filename, 'a+');
   
	for j=1:length(fpga_iterates_out)
		fprintf(fid, '%2.18f,',fpga_iterates_out(j));
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
% 	[matlab_primal_out, matlab_dual_out, matlab_aux_primal_out, matlab_iterates_out] = foo_user(project_name,state0_in, primal0_in, dual0_in, tol_iterates_in);


% 	%save matlab_primal_out_log
% 	if (TYPE_TEST==0)
% 		filename = strcat('../../ip_prototype/test/results/', project_name ,'/matlab_primal_out_log.dat');
% 	else
% 		filename = strcat('../test/results/', project_name ,'/matlab_primal_out_log.dat');
% 	end
% 	fid = fopen(filename, 'a+');
%    
% 	for j=1:length(matlab_primal_out)
% 		fprintf(fid, '%2.18f,',matlab_primal_out(j));
% 	end
% 	fprintf(fid, '\n');
% 
% 	fclose(fid);
% 
% 
% 	%save matlab_dual_out_log
% 	if (TYPE_TEST==0)
% 		filename = strcat('../../ip_prototype/test/results/', project_name ,'/matlab_dual_out_log.dat');
% 	else
% 		filename = strcat('../test/results/', project_name ,'/matlab_dual_out_log.dat');
% 	end
% 	fid = fopen(filename, 'a+');
%    
% 	for j=1:length(matlab_dual_out)
% 		fprintf(fid, '%2.18f,',matlab_dual_out(j));
% 	end
% 	fprintf(fid, '\n');
% 
% 	fclose(fid);
% 
% 
% 	%save matlab_aux_primal_out_log
% 	if (TYPE_TEST==0)
% 		filename = strcat('../../ip_prototype/test/results/', project_name ,'/matlab_aux_primal_out_log.dat');
% 	else
% 		filename = strcat('../test/results/', project_name ,'/matlab_aux_primal_out_log.dat');
% 	end
% 	fid = fopen(filename, 'a+');
%    
% 	for j=1:length(matlab_aux_primal_out)
% 		fprintf(fid, '%2.18f,',matlab_aux_primal_out(j));
% 	end
% 	fprintf(fid, '\n');
% 
% 	fclose(fid);
% 
% 
% 	%save matlab_iterates_out_log
% 	if (TYPE_TEST==0)
% 		filename = strcat('../../ip_prototype/test/results/', project_name ,'/matlab_iterates_out_log.dat');
% 	else
% 		filename = strcat('../test/results/', project_name ,'/matlab_iterates_out_log.dat');
% 	end
% 	fid = fopen(filename, 'a+');
%    
% 	for j=1:length(matlab_iterates_out)
% 		fprintf(fid, '%2.18f,',matlab_iterates_out(j));
% 	end
% 	fprintf(fid, '\n');
% 
% 	fclose(fid);

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
