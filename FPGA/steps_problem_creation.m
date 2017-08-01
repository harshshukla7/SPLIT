%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% data creation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% step 1: set the dimensions
%% step 2: change the folder in MATLAB directory
%% step 3: create a data 
%% step 4: comment so that no new data is generated
%% step 5: change the folder to PL
%% step 5a: put linear solve for FPGA
%% step 5b: put matrix vector to FPGA
%% step 6: Generate c and h file
%% step 7: change folder to SoC
%% step 7a: put linear solve for FPGA
%% step 7b: put mat-vec for suitesparse
%% step 8: rerun 
%% step 9: change folder to sw
%% step 10: select linear solve to be suitesparse
%% step 11: select matrix vector to be suitesparse
%% step 12: rerun the code


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% PL implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 1: create a folder for implementation
%% step 2: change matlab and copy paste test file
%% step 3: run the test file
%% step 4: copy paste probdata h files from PL folder
%% step 5: copy paste matrix_ops and ama .h file
%% step 6: remove Andreas code
%% step 7: copy paste foo_user till probdata
%% step 8: change h file
%% step 8a: typedef float REAL
%% step 8b: typefef float real
%% step 8c: find and replace double with float
%% step 9: copy paste data from probdata.c to foo_user.c
%% step 10: add ama.h to foo user


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% SoC implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 1: create a folder for implementation
%% step 2: change matlab and copy paste test file
%% step 3: change the dimension of input and output going to FPGA
%% step 4: copy paste fama c and h file
%% step 5: copy paste probdata c and h file
%% step 6: comment first 5 lines of test file
%% stpe 7: replace double with float in .h file 
%% step 8: add typedef float real;
%% step 9: comment declaration in probdata.c file
%% step 10: change foo_user file for nPrimal
%% step 11: soc_user add header file
%% step 12: matrix operation.c double prox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% sw implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 1: create a folder for implementation
%% step 2: change matlab and copy paste test file
%% step 3: c99 declaration
%% step 4: copy paste c data file