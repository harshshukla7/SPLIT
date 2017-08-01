function A_2D = indexes_of_nnz_elements(A)

max_nnz_elements = -1000;
indexes = {};
for i = 1:size(A, 1)
    nnz_mine = find(abs(A(i,:)) >= 1e-10);
    num_of_elements = length(nnz_mine);
    if(num_of_elements > max_nnz_elements)
       max_nnz_elements = num_of_elements;
    end
    indexes{i} = [num_of_elements nnz_mine-1];
end

% create a 2D array for export 
A_2D = zeros(size(A, 1), max_nnz_elements + 1);
for i = 1:length(indexes)
    for j = 1:length(indexes{i})
        A_2D(i,j) = indexes{i}(j);
    end
end
