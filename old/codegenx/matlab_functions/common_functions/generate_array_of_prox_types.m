function generate_array_of_prox_types(file, prox_struct, nProx)

fw = fopen(file,'a');

if nProx > 1

    fprintf(fw, 'char *prox_types[%d] = { ', nProx);
    for i = 1:nProx - 1
        fprintf(fw, '"%s"', prox_struct{i});
        fprintf(fw, ', ');
    end
    % export the last prox
    fprintf(fw, '"%s" ', prox_struct{i+1});  
    fprintf(fw, '};');
    
else
    
    fprintf(fw, 'char *prox_types[%d] = {"%s"};', nProx, prox_struct{1});
    
end

fprintf(fw, '\n');
fclose(fw);