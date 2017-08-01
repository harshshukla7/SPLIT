function generate_external_variables_matrices(file, stringM, data_type, n, m)

fw = fopen(file,'a');

fprintf(fw,'%s %s[%d][%d];', data_type ,stringM, n, m);

fprintf(fw, '\n');

fclose(fw);