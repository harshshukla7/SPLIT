%% Read text file

clear all
close all

algo = 'FPDA';
switch algo
    
    case 'FAMA'
        fileID = fopen('FAMA_2x_KKT_invert.txt','r');
        formatSpec = '%d %d %d %d %f %f %f %f %f %f';
        sizeA = [10 6];
        
        A_invert = (fscanf(fileID,formatSpec,sizeA))';
        fclose(fileID);
        
        
        fileID = fopen('FAMA_2x_KKT_ldlmat.txt','r');
        formatSpec = '%d %d %d %d %f %f %f %f %f %f';
        %sizeA = [9 Inf];
        
        A_ldlmat = (fscanf(fileID,formatSpec,sizeA))';
        fclose(fileID);
        
        fileID = fopen('FAMA_2x_KKT_ldlss.txt','r');
        formatSpec = '%d %d %d %d %f %f %f %f %f %f';
        %sizeA = [9 Inf];
        
        A_ldlss = (fscanf(fileID,formatSpec,sizeA))';
        fclose(fileID);
        
        
        fileID = fopen('FAMA_2x_KKT_ldl_lp.txt','r');
        formatSpec = '%d %d %d %d %f %f %f %f %f %f';
        %sizeA = [9 Inf];
        
        A_ldl_lp = (fscanf(fileID,formatSpec,sizeA))';
        fclose(fileID);
        
    case 'ADMM'
        
        fileID = fopen('ADMM_2x_KKT_invert.txt','r');
        formatSpec = '%d %d %d %d %f %f %f %f %f %f';
        sizeA = [10 6];
        
        A_invert = (fscanf(fileID,formatSpec,sizeA))';
        fclose(fileID);
        
        
        fileID = fopen('ADMM_2x_KKT_ldl.txt','r');
        formatSpec = '%d %d %d %d %f %f %f %f %f %f';
        %sizeA = [6 10];
        
        A_ldlmat = (fscanf(fileID,formatSpec,sizeA))';
        fclose(fileID);
        
        fileID = fopen('ADMM_2x_KKT_ldl_ss.txt','r');
        formatSpec = '%d %d %d %d %f %f %f %f %f %f';
        %sizeA = [6 10];
        
        A_ldlss = (fscanf(fileID,formatSpec,sizeA))';
        fclose(fileID);
        
        
        fileID = fopen('ADMM_2x_KKT_ldl_lp.txt','r');
        formatSpec = '%d %d %d %d %f %f %f %f %f %f';
        %sizeA = [6 10];
        
        A_ldl_lp = (fscanf(fileID,formatSpec,sizeA))';
        fclose(fileID);
        
        
    case 'FPDA'
        
        
        fileID = fopen('FPDA_2x_chol_invert.txt','r');
        formatSpec = '%d %d %d %d %f %f %f %f %f %f';
        sizeA = [10 6];
        A_invert_tp = (fscanf(fileID,formatSpec,sizeA))';
        fclose(fileID);
        A_invert = A_invert_tp;
        %A_invert = (A_invert_tp(1:32,:)+ A_invert_tp(33:64,:)+A_invert_tp(65:96,:)+A_invert_tp(97:128,:))/4;
        
        fileID = fopen('FPDA_2x_chol_ldl.txt','r');
        formatSpec = '%d %d %d %d %f %f %f %f %f %f';
        %sizeA = [9 Inf];
        A_ldlss_tp = (fscanf(fileID,formatSpec,sizeA))';
        fclose(fileID);
        A_ldlss = A_ldlss_tp;
        %A_ldlss = (A_ldlss_tp(1:32,:)+ A_ldlss_tp(33:64,:)+A_ldlss_tp(65:96,:)+A_ldlss_tp(97:128,:))/4;
        
        fileID = fopen('FPDA_2x_chol_llt.txt','r');
        formatSpec = '%d %d %d %d %f %f %f %f %f %f';
        %sizeA = [9 Inf];
        A_cholmat_tp = (fscanf(fileID,formatSpec,sizeA))';
        fclose(fileID);
        A_cholmat = A_cholmat_tp;
        %A_cholmat = (A_cholmat_tp(1:32,:)+ A_cholmat_tp(33:64,:)+A_cholmat_tp(65:96,:)+A_cholmat_tp(97:128,:))/4;
        
        fileID = fopen('FPDA_2x_chol_CLPACK_final.txt','r');
        formatSpec = '%d %d %d %d %f %f %f %f %f %f';
        %sizeA = [10 Inf];
        A_ldllp_tp = (fscanf(fileID,formatSpec,sizeA))';
        fclose(fileID);
        A_ldllp = A_ldllp_tp;
        %A_ldllp = (A_ldllp_tp(1:6,:)+ A_ldllp_tp(7:12,:))/2;
        
    case 'MAT_VEC'
        
        
end