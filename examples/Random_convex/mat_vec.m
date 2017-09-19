clear all
close all

%% Steps taken 
% 10:10:100 of the step of 5
% 100:100:1000 of the step 20

warning('this file is commented to get timing for exhaustive code generation')
band = 1;
itr = 4000; %4000;
Mat_Vec = 'exhaustive_gen'; % 'blas', 'ss', 'for_loops', sparse', 'exhaustive_gen'

%% first i = 10:10 :100 full and eye only 

% for i= 10 :10 :100
%     
%     
%     A = speye(i);
%     b = rand(i,1);
%     c = zeros(i,1);
%     
%     clear hsd;
%     hsd = splitData;
%     
%     hsd.define('itr',  itr, 'int');
%     hsd.define('bd',  0, 'int');
%     hsd.define('zerosA',nnz(~A),'int');
%     hsd.define('nonzerosA',nnz(A),'int');
%     hsd.define('totalA',numel(A),'int');
%     
%     
%     
%     a_name = sprintf('vec_b');
%     hsd.add_var(a_name, b, 'type', 'real');
%     
%     a_name = sprintf('vec_c');
%     hsd.add_var(a_name, c, 'type', 'real');
%     
%     hsd.add_function(coderFunc_times('custom_mult_A', sparse(A) , 'Astr', 'A', 'method', Mat_Vec));
%     
%     for k=1:1:itr
%         c = A*b;
%         b = A*c;
%     end
%     
%     a_name = sprintf('sol');
%     hsd.add_var(a_name, b, 'type', 'real');
%     
%     
%     hsd.write_to_file('probData');
%     
%     
%     %% Run In Terminal
%     
%     system('gcc -std=c99 -Wall -O3 -framework Accelerate -o mat test_mat_vec.c splitLoad.c probData.c splitTimer.c');
%     system('./mat');
%     
%     for j = 5 :5:i
%         
%         n = i;
%         band = j;
%         A_tp = 0.2*createbandOPT(n,band);
%         A = tril(A_tp);
%         b = rand(n,1);
%         c = zeros(n,1);
%         
%         clear hsd;
%         hsd = splitData;
%         
%         hsd.define('itr',  itr, 'int');
%         hsd.define('bd',  band, 'int');
%         hsd.define('zerosA',nnz(~A),'int');
%         hsd.define('nonzerosA',nnz(A),'int');
%         hsd.define('totalA',numel(A),'int');
%         
%         
%         
%         a_name = sprintf('vec_b');
%         hsd.add_var(a_name, b, 'type', 'real');
%         
%         a_name = sprintf('vec_c');
%         hsd.add_var(a_name, c, 'type', 'real');
%         
%         hsd.add_function(coderFunc_times('custom_mult_A', sparse(A) , 'Astr', 'A', 'method', Mat_Vec));
%         
%         
%         for k=1:1:itr
%             c = A*b;
%             b = A*c;
%         end
%         
%         
%         a_name = sprintf('sol');
%         hsd.add_var(a_name, b, 'type', 'real');
%         
%         
%         hsd.write_to_file('probData');
%         
%         
%         %% Run In Terminal
%         
%         system('gcc -std=c99 -Wall -framework Accelerate -o mat test_mat_vec.c splitLoad.c probData.c splitTimer.c');
%         system('./mat');
%         
%     end
% end


%% third i = 100:100:1000 band jump of 20 and eye as well
for i=800  %:100:1000
    
    
    A = speye(i);
    b = rand(i,1);
    c = zeros(i,1);
    
    clear hsd;
    hsd = splitData;
    
    hsd.define('itr',  itr, 'int');
    hsd.define('bd',  0, 'int');
    hsd.define('zerosA',nnz(~A),'int');
    hsd.define('nonzerosA',nnz(A),'int');
    hsd.define('totalA',numel(A),'int');
    
    
    
    a_name = sprintf('vec_b');
    hsd.add_var(a_name, b, 'type', 'real');
    
    a_name = sprintf('vec_c');
    hsd.add_var(a_name, c, 'type', 'real');
    
    
    hsd.add_function(coderFunc_times('custom_mult_A', sparse(A) , 'Astr', 'A', 'method', Mat_Vec));
    
    disp('uncooment to following to get actual timings')
%     for k=1:1:itr
%         c = A*b;
%         b = A*c;
%     end
%     
    a_name = sprintf('sol');
    hsd.add_var(a_name, b, 'type', 'real');
    
    
    hsd.write_to_file('probData');
    
    
    %% Run In Terminal
    
    disp('uncooment to following to get actual timings')
%     system('gcc -std=c99 -Wall -O3 -framework Accelerate -o mat test_mat_vec.c splitLoad.c probData.c splitTimer.c');
%     system('./mat');
    
    for j = 20 :20: i
        
        n = i;
        band = j;
        A_tp = 0.2*createbandOPT(n,band);
        A = tril(A_tp);
        b = rand(n,1);
        c = zeros(n,1);
        
        clear hsd;
        hsd = splitData;
        
        hsd.define('itr',  itr, 'int');
        hsd.define('bd',  band, 'int');
        hsd.define('zerosA',nnz(~A),'int');
        hsd.define('nonzerosA',nnz(A),'int');
        hsd.define('totalA',numel(A),'int');
        
        a_name = sprintf('vec_b');
        hsd.add_var(a_name, b, 'type', 'real');
        
        a_name = sprintf('vec_c');
        hsd.add_var(a_name, c, 'type', 'real');
        
        hsd.add_function(coderFunc_times('custom_mult_A', sparse(A) , 'Astr', 'A', 'method', Mat_Vec));
        
        disp('uncooment to following to get actual timings')
%         for k=1:1:itr
%             c = A*b;
%             b = A*c;
%         end
%         
        
        a_name = sprintf('sol');
        hsd.add_var(a_name, b, 'type', 'real');
        
        
        hsd.write_to_file('probData');
        
        
        %% Run In Terminal
        disp('uncooment to following to get actual timings')
%         system('gcc -std=c99 -Wall -framework Accelerate -o mat test_mat_vec.c splitLoad.c probData.c splitTimer.c');
%         system('./mat');
        
    end
end