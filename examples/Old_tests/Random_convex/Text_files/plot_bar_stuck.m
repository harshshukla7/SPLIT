% NumStacksPerGroup = 3;
% NumGroupsPerAxis = 6;
% NumStackElements = 4;
%
% % labels to use on tick marks for groups
% groupLabels = { 'Test'; 2; 4; 6; 8; -1; };
% stackData = rand(NumGroupsPerAxis,NumStacksPerGroup,NumStackElements);
%
% plotBarStackGroups(stackData, groupLabels);
% set(gca,'FontSize',18)
% set(gcf,'Position',[100 100 720 650])
% grid on
% set(gca,'Layer','top') % put grid lines on top of stacks

% nu = [2, 5, 10, 15, 15, 20]
% nx = [4, 10, 20, 30, 30, 40] 
% N = [4, 10, 10, 15, 20, 20] 
% (N+1)
% ADMM (N+1) nx + N(nu)
% PDA (N+1) nx 
%% 

algo = 'ADMM';

switch algo
    
    case 'ADMM'
        
        %groupLabels = { 28 ; 160 ; 320; 705; 930 ; 1240; };
        groupLabels = {48;270;540;1185;1560;2080}; % variables in KKT system
        kk = {1, 2, 3, 4, 5, 6}; % variable in the file
        
        NumStacksPerGroup = 4;
        NumGroupsPerAxis = 6;
        NumStackElements = 2;
        
        for i=1:length(kk)
            
            kkk = kk{i};
            tp (1,1:2) = [ A_ldlmat(kkk,4)/A_ldlmat(kkk,3), (A_ldlmat(kkk,5)-A_ldlmat(kkk,4))/A_ldlmat(kkk,3)];
            tp (2,1:2) = [ A_ldlss(kkk,4)/A_ldlss(kkk,3), (A_ldlss(kkk,5)-A_ldlss(kkk,4))/A_ldlss(kkk,3)];
            tp (3,1:2) = [ A_invert(kkk,4)/A_invert(kkk,3), (A_invert(kkk,5)-A_invert(kkk,4))/A_invert(kkk,3)];
            tp (4,1:2) = [ A_ldl_lp(kkk,4)/A_ldl_lp(kkk,3), (A_ldl_lp(kkk,5)-A_ldl_lp(kkk,4))/A_ldl_lp(kkk,3)];
            
            st_hs(i,:,:) =  tp;
            str ='KKT time';
            
        end
        
    case 'FAMA'
        
         %groupLabels = { 10 ; 55; 155; 310; 690 ; 920; 1220; };
        
        %kk = {17, 32, 2, 6, 11, 15, 16};
        
        %groupLabels = { 48 ; 270 ; 540; ; 930 ; 1240; };
        groupLabels = {48;270;540;1185;1560;2080};
        kk = {1, 16, 17, 22, 27,31};
        
        
        NumStacksPerGroup = 4;
        NumGroupsPerAxis = 6;
        NumStackElements = 2;
        
        for i=1:length(kk)
            
            kkk = kk{i};
            tp (1,1:2) = [ A_ldlmat(kkk,4)/A_ldlmat(kkk,3), (A_ldlmat(kkk,5)-A_ldlmat(kkk,4))/A_ldlmat(kkk,3)];
            tp (2,1:2) = [ A_ldlss(kkk,4)/A_ldlss(kkk,3), (A_ldlss(kkk,5)-A_ldlss(kkk,4))/A_ldlss(kkk,3)];
            tp (3,1:2) = [ A_invert(kkk,4)/A_invert(kkk,3), (A_invert(kkk,5)-A_invert(kkk,4))/A_invert(kkk,3)];
            tp (4,1:2) = [ A_ldl_lp(kkk,4)/A_ldl_lp(kkk,3), (A_ldl_lp(kkk,5)-A_ldl_lp(kkk,4))/A_ldl_lp(kkk,3)];
            
            st_hs(i,:,:) =  tp;
            
            str ='KKT time';
        end
        
    case 'FPDA'
         
        %groupLabels = { 10 ; 55 ; 110; 240; 315 ; 420; };
        groupLabels = {20;110;220;480;630;840};
        

        kk = {1, 16, 19, 23, 28, 32};
        
        
        NumStacksPerGroup = 4;
        NumGroupsPerAxis = 6;
        NumStackElements = 2;
        
        for i=1:length(kk)
            
            kkk = kk{i};
            tp (1,1:2) = [ A_cholmat(kkk,4)/A_cholmat(kkk,3), (A_cholmat(kkk,5)-A_cholmat(kkk,4))/A_cholmat(kkk,3)];
            tp (2,1:2) = [ A_ldlss(kkk,4)/A_ldlss(kkk,3), (A_ldlss(kkk,5)-A_ldlss(kkk,4))/A_ldlss(kkk,3)];
            tp (3,1:2) = [ A_invert(kkk,4)/A_invert(kkk,3), (A_invert(kkk,5)-A_invert(kkk,4))/A_invert(kkk,3)];
            tp (4,1:2) = [ A_ldllp(i,5)/A_ldllp(i,4),(A_ldllp(i,6)-A_ldllp(i,5))/A_ldllp(i,4)   ];
            
            st_hs(i,:,:) =  tp;
            
            str = 'Lin. Solve';
        end
end


h_tp = plotBarStackGroups(st_hs, groupLabels);
set(gca,'FontSize',18)
set(gcf,'Position',[100 100 720 650])
grid on
set(gca,'Layer','top') % put grid lines on top of stacks
xlabel('Primal Variables');
ylabel('Time (ns) per iteration');


%% Setting for Survey

hleg = legend([h_tp(1,1),h_tp(2,1),h_tp(3,1),h_tp(4,1),h_tp(1,2)],{'Custom','SuiteSparse','Invert','LAPACK','Remainders'});
hcg = title('PDA');
set(hcg, 'FontSize', 16);
ylabcg = ylabel('Time (ns) per iteration');
set(ylabcg, 'FontSize', 16);
xlabcg = xlabel('No. of variables');
set(xlabcg, 'FontSize', 16);
set(hleg,'FontSize', 14);
set(gca, 'FontName', 'Helvetica')
set(ylabcg, 'FontName', 'Helvetica')
set(xlabcg, 'FontName', 'Helvetica')
set(hleg, 'FontName', 'Helvetica')
grid off



%%
%legend(str,'Rest Operations')
% if strcmp(algo,'FPDA')
%  legend([h_tp(1,1),h_tp(2,1),h_tp(3,1),h_tp(1,2)],{'Custom','SuiteSparse','Invert','LAPACK','Remainders'})   
% else
% legend([h_tp(1,1),h_tp(2,1),h_tp(3,1),h_tp(4,1),h_tp(1,2)],{'Custom','SuiteSparse','Invert','Lapack','Remainders'})
% end
% legend([h_tp(1,1),h_tp(2,1),h_tp(3,1),h_tp(4,1),h_tp(1,2)],{'hsf','isf','tsf','bsf','gsf'})

%%  Plot 3D 

% mat = A_cholmat ;
% 
% nxu_vec = zeros(32,1);
% N_vec = zeros(32,1);
% value_vec = zeros(32,1);
% i = 1;
% for nxu = 2:1:5
%     for N = 2:1:5
%         indx = find(mat(:,1) == (2*nxu*N) + nxu);
%         value = mat(indx(1),4)/mat(indx(1),3);
%         nxu_vec(i,1) = 2*nxu;
%         N_vec(i,1) = N;
%         value_vec(i,1) = value;
%         i = i+1;
%     end 
% end
%         
% for nxu = 5:5:20
%     for N = 5:5:20
%         indx = find(mat(:,1) == (2*nxu*N) + nxu);
%         value = mat(indx(1),4)/mat(indx(1),3);
%         nxu_vec(i,1) = 2*nxu;
%         N_vec(i,1) = N;
%         value_vec(i,1) = value;
%         i = i+1;
%     end 
% end
%         
        
        
