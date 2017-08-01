function generate_matrix(file, M, stringM, data_type, n, m)

fw = fopen(file,'a');

fprintf(fw,'%s %s[%d][%d] = {', data_type ,stringM, n, m);
for i = 1:n
    fprintf(fw,' {');
    for j = 1:m
        if(j == m)
            fprintf(fw,'%d ', M(i,j));
        else
            fprintf(fw,'%d, ',M(i,j)); 
        end
    end
    if(i == n)
        fprintf(fw,'} ');
    else
        fprintf(fw,'}, ');
    end
end
fprintf(fw,'}; \n');

fclose(fw);