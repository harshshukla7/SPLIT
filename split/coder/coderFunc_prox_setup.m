classdef coderFunc_prox_setup < coderFunc
    
    methods
        function f = coderFunc_prox_setup(dat,setting)
            % Create a function that evaluates the prox functions for the data
            %
            % y = prox(x)
            %
            
            f = f@coderFunc('void custom_prox_setup(double y[nDual], double x[nDual])');
            
            
            
        end
    end
end