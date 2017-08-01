classdef coderFunc_mldivide_chol < coderFunc
    
    methods
        function f = coderFunc_mldivide(func_name, Q, A, varargin)
            %
            % func_mldivide(dat, func_name, A, param1, val1, etc)
            %
            % Generate code to solve a system of equations
            %
            % y = A\x
            %
            % Function prototype:
            %   void func_name(y[n], x[n])
            %
            % Arguments:
            %  func_name  : string
            %  A          : Matrix data
            %
            % Optional:
            %  Astr       : Name to use when storing the matrix [otherwise random]
            %  method     : One of 'ldl', 'invert'
            %               If not specified, method will be chosen based on the
            %               structure of A
            %  desc       : Text description of the function
            %  zero_tol   : Value below which a number is considered zero [1e-10]
            %               Used to sparsify the invert method
            %
            
            p = inputParser;
            addRequired(p, 'func_name',  @ischar);
            addRequired(p, 'A',          @isnumeric);
            addRequired(p, 'Q',          @isnumeric);
            
            addParameter(p, 'Astr',       uniqueVarName, @ischar);
            addParameter(p, 'method', 'auto', ...
                @(x) any(validatestring(x,{'choles'})));
            addParameter(p, 'desc', '', @ischar);
            addParameter(p, 'zero_tol', 1e-10, @isnumeric);
            
            % Parse and copy the results to the workspace
            parse(p, func_name, A, varargin{:});
            cellfun(@(q) evalin('caller',[q ' = p.Results.' q ';']), fieldnames(p.Results))
            
            % Write the function prototype
            %             [m,n] = size(A);
            %             assert(m==n, 'func_mldivide currently only solves square matrices')
            f = f@coderFunc('void %s(REAL y[%i], REAL x[%i])', func_name, n, n);
            
            
            
            
            method = 'choles';
            
            f.desc = method;
            
            
            n_KKT = size(Q,1) + size(A,2);
            
            
            f.add_func(coderFunc_prefactor('custom_compute_prefactor',A,'none'));
            
            % step 1: add variables, K21K11^(-1), K21, K11^(-1) cholesky
            % K11
            %%%%%% matrix A
            rind_name = sprintf('A_chol_rind');
            cind_name = sprintf('A_chol_cind');
            val_name = sprintf('A_chol_val');
            
            [rind_L_ret, cind_L_ret, val_L_ret] = L_mat_vec_split(A);
            len_val = length(val_L_ret);
            out_size_A = size(A,1);
            
            f.add_var(rind_name, rind_L_ret, 'type', 'int');
            f.add_var(cind_name, cind_L_ret, 'type', 'int');
            f.add_var(val_name, val_L_ret, 'type', 'real');
            
            % step 2: set another vector entries to zeros
            f.pl('int i, j;');
            f.pl('float x_int[%d]', n_KKT);
            
            % step 3: mat-vec
            
            
            % step 4: vector addition
            % step 5: lower cholesky solve
            % step 6: upper cholesky solve
            % step 7: mat-vec
            
            
            
            
        end
    end
end


