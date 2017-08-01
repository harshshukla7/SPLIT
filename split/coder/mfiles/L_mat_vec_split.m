function [rind_L_ret, cind_L_ret, val_L_ret] = L_mat_vec_split(L)

[rows_L, columns_L] = size(L);

[rind_L, cind_L, val_L] = find(L);




max_nnz_rows = 0;
min_nnz_rows = rows_L;

% %  first we find maximum and minimum number of non zeros per column in
% entire matrix

for i = 1:columns_L
    
    tmp = nnz(L(:,i));
    
    if tmp>max_nnz_rows
        
        max_nnz_rows = tmp;
    end
    
    if tmp<min_nnz_rows
        
        min_nnz_rows = tmp;
    end
    
end

non_zero_columns = 0;
for i=1:columns_L
    
    
    tmp = nnz(L(:,i));
    
    if (tmp == max_nnz_rows) == 1
        
        non_zero_columns = non_zero_columns + 1;
        
        
    elseif (tmp == 0) == 1
        
    else
        
        non_zero_columns = non_zero_columns + 1;
        need_fill = max_nnz_rows - tmp;
        
        % step 1 find index of nonzeros
        % step 2 fill zeros in values and row and column index
        
        [tmp_row, ~, ~] = find(L(:,i));
        
        for j = 1:need_fill
            
            for k = 1:rows_L
                if isempty(find(tmp_row == k)) == 1
                    
                    
                    rind_L = [rind_L; k];
                    cind_L = [cind_L; i];
                    val_L = [val_L; 0];
                    
                    break;
                    
                    
                end
            end
            
            
        end
        
    end
    
    
    
end

% %  if we do not have the maximum and minimum number of nonzeros the same value then we
% need to take care of this smartly and current "smartly" is to give an
% error!!!!! :P

% if max_nnz_rows ~= min_nnz_rows
%
%     error('currently SPLIT does not support varying number of non zeros per columns ')
%
% end

total_entry = max_nnz_rows*non_zero_columns;
rind_L_ret = zeros(total_entry,1);
cind_L_ret = zeros(total_entry,1);
val_L_ret = zeros(total_entry,1);


start_ind = -max_nnz_rows + 1;
for j = 1:columns_L
    
    [i_col_ind_tmp, act_ind] = find(cind_L==j); % find all entries of ith row
    
    if isempty(i_col_ind_tmp) == 1
        
    else
        %       i_row_ind_tmp = rind_L(act_ind); % find all columns of ith row
        
        start_ind = start_ind + max_nnz_rows ;
        end_ind   = start_ind + max_nnz_rows - 1;
        
        true_indx = i_col_ind_tmp;
        rind_L_ret(start_ind:end_ind,1) = rind_L(true_indx);
        cind_L_ret(start_ind:end_ind,1) = cind_L(true_indx);
        val_L_ret(start_ind:end_ind,1) = val_L(true_indx);
    end
    
end


end