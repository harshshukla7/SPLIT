classdef coderFunc_dual_update < coderFunc
    
    methods
        function f = coderFunc_dual_update(dat,setting)
            % Create a function that evaluates the prox functions for the data
            %
            % y = prox(x)
            %
            
            f = f@coderFunc('void custom_dual_update(double y[nDual], double x[nDual])');
            
            
            
        end
    end
end