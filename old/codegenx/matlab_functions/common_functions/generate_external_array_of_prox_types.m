function generate_external_array_of_prox_types(file)

fw = fopen(file,'a');

fprintf(fw, 'extern char prox_types[][40]; ');

fprintf(fw, '\n');

fclose(fw);