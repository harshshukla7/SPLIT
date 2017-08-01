function generate_vector(file, M, stringM, data_type, n)

fw = fopen(file,'a');

fprintf(fw,'%s %s[%d] = {',data_type, stringM, n);
for i = 1:n
    if(i == n)
       fprintf(fw,'%d ' , M(i));
    else
       fprintf(fw,'%d, ', M(i)); 
    end
end
fprintf(fw,'}; \n');

% IMPORTANT: never forget to close the corresponding file
% descriptor, because then I get error messagees: std::exception !!!!!
fclose(fw);