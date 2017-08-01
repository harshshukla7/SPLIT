function generate_external_variables(file, stringM, data_type, n)

fw = fopen(file,'a');

fprintf(fw,'%s %s[%d]; ',data_type, stringM, n);

fprintf(fw, '\n');

% IMPORTANT: never forget to close the corresponding file
% descriptor, because then I get error messagees: std::exception !!!!!
fclose(fw);