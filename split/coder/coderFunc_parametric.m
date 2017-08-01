classdef coderFunc_parametric < coderFunc
    
    methods
        function func = coderFunc_parametric(fname, dat, varargin)
            %
            % f = coderFunc_parametric(fname, dat)
            %
            % Generates a function to compute the parametric values
            %   f = pF*par + f_const
            %   l = pL*par + l_const
            %   b = pB*par + b_const
            %
            %  void fname(REAL l[xxx], REAL f[xxx], REAL b[xxx], REAL par[xxx])
            %
            
            
            %             addParameter(p, 'method', 'auto', ...
            %                 @(x) any(validatestring(x,{'auto','for_loops','sparse','FPGA_matvec', 'FPGA_matvec_dense','blas','ss','exhaustive_gen'})));
            
             p = inputParser;
            addParameter(p, 'method', 'auto', ...
                @(x) any(validatestring(x,{'auto', 'auto_FPGA', 'for_loops','sparse','FPGA_matvec', 'FPGA_matvec_dense', 'FPGA_exhaustive_gen','blas','ss','exhaustive_gen'})));
            
            % Parse and copy the results to the workspace
            parse(p, varargin{:});
            cellfun(@(q) evalin('caller',[q ' = p.Results.' q ';']), fieldnames(p.Results))
            
            
            func = func@coderFunc('void %s(REAL l[%i], REAL f[%i], REAL b[%i], REAL par[%i])',...
                fname, length(dat.l), length(dat.f), length(dat.b), size(dat.pL,2));
            
            vars = {{'pL', 'l'}, {'pF', 'f'}, {'pB', 'b'}};
            for i = 1:length(vars)
                pVar = vars{i}{1};
                cVar = vars{i}{2};
                
                if any(abs(dat.(pVar)(:))) > 0
                    % We need to compute the product, since pL is non-zero
                    func.add_func(coderFunc_times(['compute_parametric_' pVar], dat.(pVar), 'bConst', dat.(cVar), 'method', method));
                    func.pl(['compute_parametric_' pVar '(' cVar ', par);']);
                    
                    % Allocate memory for the variable b
                    func.add_var(cVar, zeros(length(dat.(cVar)),1), 'type', 'real');
                    
                else
                    % pL is zero, so l is a constant - just store it and do nothing
                    func.add_var(cVar, dat.(cVar), 'storage_method', 'dense');
                end
            end
        end
    end
end
