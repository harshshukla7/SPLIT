clear all
close all



%%

fileID = fopen('MAT_VEC_text_blas.txt','r');
formatSpec = '%d %d %d %d %d %d %f %f %f %f';
sizeA = [10 Inf];
A_blas = (fscanf(fileID,formatSpec,sizeA))';
fclose(fileID);


fileID = fopen('MAT_VEC_text_for_loops.txt','r');
formatSpec = '%d %d %d %d %d %d %f %f %f %f';
sizeA = [10 Inf];
A_for_loops = (fscanf(fileID,formatSpec,sizeA))';
fclose(fileID);

fileID = fopen('MAT_VEC_text_ss.txt','r');
formatSpec = '%d %d %d %d %d %d %f %f %f %f';
sizeA = [10 Inf];
A_ss = (fscanf(fileID,formatSpec,sizeA))';
fclose(fileID);


fileID = fopen('MAT_VEC_text_exchustive.txt','r');
formatSpec = '%d %d %d %d %d %d %f %f %f %f';
sizeA = [10 Inf];
A_exhaustive = (fscanf(fileID,formatSpec,sizeA))';
fclose(fileID);

%% 

full_index = find ((A_blas(1:354,1)./A_blas(1:354,2)==1));
full_index = full_index([1:10,12:end]);
half_index = find(A_blas(1:354,1)./A_blas(1:354,2)==2);
half_index = [half_index;145;187;249;331];
half_index = sort(half_index);
diag_index = find (A_blas(1:354,2)==0);
diag_index = diag_index([1:10,12:end]);
size_index = A_blas(diag_index,1);
%size_index = size_index([1:10,12:end]);


%% color definition
ss_c=[230	69	69	]/255; %% SuiteSparse
ct_c=[37 178 207]/255; %% Custom
it_c=[77 200 161]/255; %% Invert
fancyyellow=[232 205 35]/255;
orange=[255 175 42]/255;
darkorange=[200 100 50]/255;
fancycyan=[0.2 0.7 0.9];
bs_c=[220 100 50]/255; %% BLAS



%% Full Index

figure()
h1 = semilogy(size_index,A_blas(full_index,7),'--','color',bs_c,'LineWidth', 2.5);
hold on
h2 = semilogy(size_index,A_exhaustive(full_index,7),'--','color',ct_c,'LineWidth', 2.5);
hold on
h3 = semilogy(size_index,A_for_loops(full_index,7),'--','color',it_c,'LineWidth', 2.5);
hold on
h4 = semilogy(size_index,A_ss(full_index,7),'--','color',ss_c,'LineWidth', 2.5);

%% half

%figure()
h5 = semilogy(size_index,A_blas(half_index,7),'-*','color',bs_c,'LineWidth', 2.5);
hold on
semilogy(size_index,A_exhaustive(half_index,7),'-*','color',ct_c,'LineWidth', 2.5)
hold on
semilogy(size_index,A_for_loops(half_index,7),'-*','color',it_c,'LineWidth', 2.5)
hold on
semilogy(size_index,A_ss(half_index,7),'-*','color',ss_c,'LineWidth', 2.5)
hold on

%% diag

%figure()
h6 = semilogy(size_index,A_blas(diag_index,7),'-+','color',bs_c,'LineWidth', 2.5);
hold on
semilogy(size_index,A_exhaustive(diag_index,7),'-+','color',ct_c,'LineWidth', 2.5)
hold on
semilogy(size_index,A_for_loops(diag_index,7),'-+','color',it_c,'LineWidth', 2.5)
hold on 
semilogy(size_index,A_ss(diag_index,7),'-+','color',ss_c,'LineWidth', 2.5)
hold on 

%%

h6 = plot(nan,nan,'color',bs_c);
hold on 
h7 = plot(nan,nan,'color',ct_c);
hold on 
h8 = plot(nan,nan,'color',it_c);
hold on 
h9 = plot(nan,nan,'color',ss_c);
hold on 
h10 = plot(nan,nan,'-','color','k');
hold on
h11 = plot(nan,nan,'*','color','k');
hold on 
h12 = plot(nan,nan,'`+','color','k');
hold on 

%%
legend([h7,h9,h8,h6,h10,h11,h12],{'Custom','SuiteSparse','Invert','BLAS','Full','Half','Eye'})
