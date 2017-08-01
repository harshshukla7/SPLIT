function generate_external_array_of_norm_types(file)

fw = fopen(file,'a');

fprintf(fw, 'extern char norm_types[][5]; ');

fprintf(fw, '\n');

fclose(fw);