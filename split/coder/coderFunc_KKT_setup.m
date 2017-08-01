classdef coderFunc_KKT_setup < coderFunc
    
    methods
        function f = coderFunc_KKT_setup(dat,setting)
            % Create a function that evaluates the prox functions for the data
            %
            % y = prox(x)
            %
            
            f = f@coderFunc('void custom_KKT_setup(double y[nDual], double x[nDual])');
            
            
            
        end
    end
end