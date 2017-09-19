clear all
close all

%% Steps taken 
% 10:10:100 of the step of 5
% 100:100:1000 of the step 20

warning('this file is commented to get timing for exhaustive code generation')
band = 1;
itr = 4000; %4000;
Mat_Vec = 'exhaustive_gen'; % 'blas', 'ss', 'for_loops', sparse', 'exhaustive_gen'

timing = [];
%% first i = 10:10 :100 full and eye only 

for i= 10 :10 :100
    i
    
    for j = i/2
        
        n = i;
        band = j;
        A_tp = 0.2*createbandOPT(n,band);
        A = tril(A_tp);
        b = rand(n,1);
        c = zeros(n,1);
        
        clear hsd;
        hsd = splitData;
        
        tic
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
        
        
        a_name = sprintf('sol');
        hsd.add_var(a_name, b, 'type', 'real');
        
        
        hsd.write_to_file('probData');
        el = toc;
        timing = [timing; i,el];
        
      
        
    end
end


%% third i = 100:100:1000 band jump of 20 and eye as well
for i=100  :100: 1000
    i
    
    for j = i/2
        
        n = i;
        band = j;
        A_tp = 0.2*createbandOPT(n,band);
        A = tril(A_tp);
        b = rand(n,1);
        c = zeros(n,1);
        
        clear hsd;
        hsd = splitData;
        
        tic
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
        
        a_name = sprintf('sol');
        hsd.add_var(a_name, b, 'type', 'real');
        
        
        hsd.write_to_file('probData');
        
        
        el = toc;
       timing = [timing; i,el];
    end
end