function generate_array_of_norm_types(file, prox_norm, nProx)

fw = fopen(file,'a');

if nProx > 1

    fprintf(fw, 'char *norm_types[%d] = { ', nProx);
    for i = 1:nProx - 1
        fprintf(fw, '"%s"', prox_norm{i});
        fprintf(fw, ', ');
    end
    % export the last prox
    fprintf(fw, '"%s" ', prox_norm{i+1});  
    fprintf(fw, '};');
    
else
    
    fprintf(fw, 'char *norm_types[%d] = {"%s"};', nProx, prox_norm{1});
    
end

fprintf(fw, '\n');
fclose(fw);