function [] = split_MV_sparse(H, file_name, reset_enable)
    % this function generates code for sparse matrix vector multiplication
    % y_out = (y_out) + H*x_in 
    % H is a sparse matrix with at least one nonzero element
    % matrix_name must be a string, it will be used for generated files and
    % functions namings
    % if reset_enable == 1 => y_out = H*x_in, else => y_out = y_out + H*x_in
 
   
    data_t = 'data_t_primal_out';

    % identify the sparsity pattern
    [row, col, val] = find(H);
    coo.row = row';
    coo.col = col';
    coo.val = val';

    % perform scheduling
    [coo_scheduled, d_max ] = schedule_mat_vec( coo );
    
    % test equivalence before and after scheduling
    %max(max(abs(sort([coo.row; coo.col; coo.val],2) - sort([coo_scheduled.row; coo_scheduled.col; coo_scheduled.val],2))));
    
    % generate header file
    tp_name = strcat('user_', file_name, '_sparse_mv_mult.h');
    fileID = fopen(sprintf(tp_name),'w');
    fprintf(fileID,'#include "foo_data.h" \n');
    fprintf(fileID,'#define %s_SIZE_row %d\n', file_name, max(row));
    fprintf(fileID,'#define %s_SIZE_col %d\n', file_name, max(col));
    fprintf(fileID,'#define %s_NNZ %d\n\n', file_name, max(size(val)));
    fprintf(fileID,strcat('void',32, file_name, '_sparse_mv_mult', '(',data_t,' y_out[',file_name,'_SIZE_row],',data_t,' x_in[',file_name,'_SIZE_col]);\n'));
    fclose(fileID);
    
    % generate cpp file
    tp_name = strcat('user_', file_name, '_sparse_mv_mult.cpp');
    fileID = fopen( sprintf(tp_name), 'w');
    fprintf(fileID,strcat('#include',32,'"user_',file_name, '_sparse_mv_mult.h"\n'));
    fprintf(fileID,'\n');
    fprintf(fileID,strcat('void',32, file_name, '_sparse_mv_mult', '(',data_t,' y_out[',file_name,'_SIZE_row],',data_t,' x_in[',file_name,'_SIZE_col])\n'));
    fprintf(fileID,'{\n');
    fprintf(fileID,strcat('\tint i;\n'));
    fprintf(fileID,strcat('\t', 'int', ' col[',file_name,'_NNZ] = {',sprintf('%d,' , reshape((coo_scheduled.col-1).',[],1)),'};\n'));
    fprintf(fileID,strcat('\t', 'int', ' row[',file_name,'_NNZ] = {',sprintf('%d,' , reshape((coo_scheduled.row-1).',[],1)),'};\n'));
    fprintf(fileID,strcat('\t', data_t, ' val[',file_name,'_NNZ] = {',sprintf('%d,' , reshape(coo_scheduled.val.',[],1)),'};\n\n'));
    
    if reset_enable
        fprintf(fileID,'\t// reset output\n');
        fprintf(fileID,strcat('\tfor(i = 0; i <',32,file_name,'_SIZE_row; i++)\n'));
        fprintf(fileID,'\t{\n');
        fprintf(fileID,strcat('\t\t','#pragma HLS PIPELINE\n'));
        fprintf(fileID,'\t\ty_out[i] = 0;\n');
        fprintf(fileID,'\t}\n\n');
    end
      
    fprintf(fileID,'\t// perform mat vec mult\n');
    fprintf(fileID,strcat('\tfor(i = 0; i < ',32,file_name,'_NNZ; i++)\n'));
    fprintf(fileID,'\t{\n');
    fprintf(fileID,strcat('\t\t','#pragma HLS DEPENDENCE variable=y_out inter distance=',num2str(d_max),' true\n'));
    fprintf(fileID,strcat('\t\t','#pragma HLS PIPELINE II=2 \n'));
    fprintf(fileID,'\t\ty_out[row[i]] += val[i]*x_in[col[i]];\n');
    fprintf(fileID,'\t}\n');
    
    fprintf(fileID,'}\n');
    
    
    
    

end


function [ coo_scheduled, d_max] = schedule_mat_vec( coo )

% determine matrix dimensions (it is assumed that every row and column has 
% at least one nonzero entry)
n_rows = max(coo.row);
n_col  = max(coo.col);

% number of nonzero entries
N = max(size(coo.val));

% sort entries according to the number of nonzeros per row
nnz_per_row = zeros(1, n_rows);
for i = 1:n_rows
    nnz_per_row(i) = sum(coo.row == i);
end
n = max(nnz_per_row); % max number of elements per row
if n > 1
    m = sum(nnz_per_row == n); % number of rows with max number of elements
    [~,I] = sort(nnz_per_row,'descend');
    coo_sorted.row = zeros(1,N);
    coo_sorted.col = zeros(1,N);
    coo_sorted.val = zeros(1,N);
    pointer = 1;
    for i = I
        current_row_indeces = [coo.row == i];
        offset = sum(current_row_indeces);
        coo_sorted.row(pointer:pointer+offset-1) = coo.row(current_row_indeces);
        coo_sorted.col(pointer:pointer+offset-1) = coo.col(current_row_indeces);
        coo_sorted.val(pointer:pointer+offset-1) = coo.val(current_row_indeces);
        pointer = pointer + offset ;
    end

    d_max = floor((N-m)/(n-1));

    % create an array of slot structures
    time_slots = [];
    slot_counter = 1;
    for i = 1:rem((N-m),(n-1)) % first type of slots
        slot.full = 0;
        slot.size = d_max + 1;
        slot.beginning = 1;
        slot.schedule = zeros(1, d_max+1);
        time_slots = [time_slots slot];
        slot_counter = slot_counter + 1;
    end
    for i = rem((N-m),(n-1))+1:n-1 % second type of slots
        slot.full = 0;
        slot.size = d_max;
        slot.beginning = 1;
        slot.schedule = zeros(1, d_max);
        time_slots = [time_slots slot];
        slot_counter = slot_counter + 1;
    end
    slot.full = 0; % third type of slots
    slot.size = m;
    slot.beginning = 1;
    slot.schedule = zeros(1, m);
    time_slots = [time_slots slot];


    % perform scheduling
    flag = 1;
    j = 1;
    while(flag)
        for i=1:n
            if(~time_slots(i).full) % if the currecnt slot is not full 
                time_slots(i).schedule(time_slots(i).beginning) = j; % schedule the current entry
                j = j + 1;
                time_slots(i).beginning = time_slots(i).beginning + 1; % increment beginning
                if(time_slots(i).beginning > time_slots(i).size) % check if the current slot is full
                    time_slots(i).full = 1;
                end
                if j > N
                    flag = 0;
                    break;
                end
            end
        end
    end


    coo_scheduled.row = coo_sorted.row([time_slots.schedule]);
    coo_scheduled.col = coo_sorted.col([time_slots.schedule]);
    coo_scheduled.val = coo_sorted.val([time_slots.schedule]);

    % check if the solution is correct
    for i=1:n_rows
        indeces_same_row = find(coo_scheduled.row == (i) );
        indeces_same_row = sort(indeces_same_row);
        for j = 1:(max(size(indeces_same_row))-1)
            if( abs(indeces_same_row(j) - indeces_same_row(j+1)) <  d_max)
                error('Error: scheduling algorithm failed!');
            end
        end
    end
    
else
    % no need to schedule
    coo_scheduled = coo;
    d_max = 10000;
end

end
