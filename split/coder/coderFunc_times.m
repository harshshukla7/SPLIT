classdef coderFunc_times < coderFunc
    
    methods
        function f = coderFunc_times(func_name, A, varargin)
            %
            % f = coderFunc_times(func_name, A, param1, val1, etc)
            %
            % Generate code to multiply a constant matrix and a vector together
            %
            % y = A*x + b
            %
            % Function prototype:
            %   void func_name(y[n], x[n])
            % or
            %   void func_name(y[n], x[n], b[n])
            %
            % Arguments:
            %  func_name  : string
            %  A          : Matrix data
            %
            % Optional:
            %  Astr       : Name to use when storing the matrix [otherwise random]
            %  bConst     : Specify if b is a constant [double]
            %  bVar       : Specify if b is a variable [logical]
            %  method     : One of 'for_loops', 'sparse', 'blas', 'exhaustive_gen'
            %               If not specified, method will be chosen based on the
            %               structure of A
            %  desc       : Text description of the function
            %  x_name     : Name used for x [x]
            %  y_name     : Name used for y [y]
            %  method     : 'forloops', 'sparse', 'blas', 'ss', exhaustive_gen'
            
            p = inputParser;
            addRequired(p, 'func_name',  @ischar);
            addRequired(p, 'A',          @isnumeric);
            
            addParameter(p, 'Astr',       uniqueVarName, @ischar);
            addParameter(p, 'bConst',  [],  @isnumeric);
            addParameter(p, 'bVar',    false,  @islogical);
            addParameter(p, 'method', 'auto', ...
                @(x) any(validatestring(x,{'auto', 'auto_FPGA', 'FPGA_diag', 'for_loops','sparse','FPGA_matvec', 'FPGA_matvec_sparse', 'FPGA_matvec_dense', 'FPGA_exhaustive_gen','blas','ss','exhaustive_gen'})));
            addParameter(p, 'desc', '', @ischar);
            addParameter(p, 'x_name', 'x', @ischar);
            addParameter(p, 'y_name', 'y', @ischar);
            addParameter(p, 'loop_reset', 1,  @isnumeric);
            
            
            % Parse and copy the results to the workspace
            parse(p, func_name, A, varargin{:});
            cellfun(@(q) evalin('caller',[q ' = p.Results.' q ';']), fieldnames(p.Results))
            
            % Write the function prototype
            [m,n] = size(A);
            proto = sprintf('void %s(REAL %s[%i], REAL %s[%i]', func_name, y_name, m, x_name, n);
            bStr = [Astr '_b'];
            if bVar
                proto = sprintf('%s, %s[%i])', proto, bStr, m);
            else
                proto = sprintf('%s)', proto);
            end
            f = f@coderFunc(proto);
            
            if strcmpi(method, 'auto')
                % Try and guess the right method
                if (m*n<40000)
                    if (nnz(A)/(m*n)<0.5)
                        method = 'exhaustive_gen';
                        
                    else
                        method = 'blas';
                    end
                else
                    method = 'blas';
                end
            end
            
            f.desc = desc;
            
            
            
            % Check if it is for FPGA
            FPGA_gen = 0;
            if (strcmp(method, 'FPGA_diag') ||  strcmp(method, 'FPGA_matvec_sparse') || strcmp(method, 'FPGA_matvec') || strcmp(method, 'FPGA_matvec_dense') || strcmp(method, 'auto_FPGA') || strcmp(method, 'FPGA_exhaustive_gen'))
                
                FPGA_gen = 1;
            end
            
            
            if strcmpi(method, 'auto_FPGA')
                % Try and guess the right method
                
                if isdiag(A) == 1
                    
                    method = 'FPGA_matvec';
                    
                    
                elseif (m*n<40000)
                    if (nnz(A)/(m*n)<0.5)
                        method = 'FPGA_exhaustive_gen';
                        
                    else
                        method = 'FPGA_matvec';
                    end
                else
                    method = 'FPGA_matvec';
                end
            end
            
            
            % Set the output equal to b
            if ~isempty(bConst)
                bConst = full(bConst);
                f.add_var(bStr, bConst);
            end
            
            if loop_reset == 1
                if ~isempty(bConst) || bVar
                    
                    if (FPGA_gen == 1)
                        
                        f.pl('int i;');
                        val_name = sprintf('%s_val',Astr);
                        f.pl('loop_%s_reset: for (i = 0; i < %d; i++)', val_name, size(A,1) );
                        f.pl('{');
                        f.pl('#pragma    HLS PIPELINE');
                        f.pl('%s[i] = %s[i];',y_name, bStr);
                        f.pl('}');
                        
                    elseif  (strcmp(method, 'FPGA_matvec_sparse') == 1)
                        
                        error('sparse matrix vector with y = Ax is supported. Do not use for y = Ax + b')
                        
                    else
                        f.pl('memcpy(%s, %s, %i*sizeof(REAL));', y_name, bStr, size(A,1));
                        
                    end
                    
                else
                    
                    if (strcmp(method, 'FPGA_matvec_sparse') == 1)
                        
                    elseif((FPGA_gen == 1) && (strcmp(method, 'FPGA_matvec_sparse') == 0))
                        
                        f.pl('int i;');
                        val_name = sprintf('%s_val',Astr);
                        f.pl('loop_%s_reset: for (i = 0; i < %d; i++)', val_name, size(A,1) );
                        f.pl('{');
                        f.pl('#pragma    HLS PIPELINE');
                        f.pl('%s[i] = 0.0;',y_name);
                        f.pl('}');
                        
                    else
                        f.pl('memset(%s, 0, %i*sizeof(REAL));', y_name, size(A,1));
                        
                    end
                end
            end
            
            
            f.desc = method;
            
            % Write the code to actually solve the problem
            switch method
                %%
                case 'for_loops'
                    % Multiply static matrix and vector using dense for-loops
                    f.add_var(Astr, full(A)');
                    
                    f.pl('for(int r=0; r<%i; r++) ', size(A,1))
                    f.pl('  for(int c=0; c<%i; c++) %s[r] += %s[r*%i+c]*%s[c];', size(A,2), y_name, Astr, size(A,2), x_name)
                    
                    %%
                case 'sparse'
                    % Write data in compressed format
                    f.add_var(Astr, A, 'type', 'real', 'storage_method', 'sparse');
                    
                    f.pl('int i=0;')
                    f.pl('for(int r=0; r<%i; r++) ', size(A,1))
                    f.pl('  for(int ic=0; ic<%s_nz_per_row[r]; ic++) {', Astr)
                    f.pl('    %s[r] += %s_dat[i]*%s[%s_ind[i]];', y_name, Astr, x_name, Astr)
                    f.pl('    i++;')
                    f.pl('  }')
                    
                    %%
                case 'blas'
                    % Multiply static matrix and vector using BLAS
                    f.add_var(Astr, full(A)');
                    
                    warning('realType must be double to use BLAS')
                    f.pl('cblas_dgemv(CblasRowMajor, CblasNoTrans, %i, %i, 1.0, (double*)%s, %i, %s, 1, 1.0, %s, 1);',...
                        size(A,1), size(A,2), Astr, size(A,2), x_name, y_name);
                    
                    %%
                case 'ss'
                    
                    warning('Current version for sparse matrix vector multiplication does not support matrix size larger than int range')
                    n_name = sprintf('n_%s_ss',Astr);
                    Lp_ss_name = sprintf('%s_p_ss',Astr);
                    Li_ss_name = sprintf('%s_i_ss',Astr);
                    Lx_ss_name = sprintf('%s_x_ss',Astr);
                    
                    %%%%% Create in CSR format
                    
                    A = sparse(A);
                    [Lp, Li, Lx] = sparse_to_csr(A'); %%% Suitesparse has column format
                    
                    Lp = Lp-1; %%% to bring into c format
                    Li = Li-1; %%% to bring into c format
                    
                    %%%%%% CSR format created
                    
                    
                    %%%% Save variable in binary file
                    f.add_var(Lp_ss_name,Lp,'type','int');
                    f.add_var(Li_ss_name,Li,'type','int');
                    f.add_var(Lx_ss_name,Lx,'type','real');
                    f.add_var(n_name,size(A,2),'type','int');
                    
                    %%%%
                    
                    f.pl('int j,p;');
                    f.pl('for (j = 0 ; j < %s[0] ; j++)', n_name);
                    f.pl('{');
                    f.pl('for (p = %s[j] ; p < %s[j+1] ; p++)',Lp_ss_name, Lp_ss_name);
                    f.pl('{');
                    f.pl('%s[%s[p]] += %s[p] * %s[j];',y_name, Li_ss_name, Lx_ss_name, x_name);
                    f.pl('}');
                    f.pl('}');
                    
                    
                case 'FPGA_matvec'
                    
                    rind_name = sprintf('%s_rind',Astr);
                    cind_name = sprintf('%s_cind',Astr);
                    val_name = sprintf('%s_val',Astr);
                    
                    [rind_L_ret, cind_L_ret, val_L_ret] = L_mat_vec_split(A);
                    len_val = length(val_L_ret);
                    out_size = size(A,1);
                    
                    rind_L_ret = rind_L_ret - 1;
                    cind_L_ret = cind_L_ret -1;
                    f.add_var(rind_name, rind_L_ret, 'type', 'int');
                    f.add_var(cind_name, cind_L_ret, 'type', 'int');
                    f.add_var(val_name, val_L_ret, 'type', 'real');
                    
                    f.pl('int row_ind, col_ind;');
                    %                     f.pl('loop_%s_reset: for (i = 0; i < %d; i++)', val_name, out_size);
                    %                     f.pl('{');
                    %                     f.pl('#pragma    HLS PIPELINE');
                    %                     f.pl('%s[i] = 0.0;',y_name);
                    %                     f.pl('}');
                    
                    f.pl('loop_%s: for (i = 0; i < %d; i++)', val_name, len_val);
                    f.pl('{');
                    f.pl('#pragma    HLS PIPELINE');
                    f.pl('row_ind = %s[i];', rind_name);
                    f.pl('col_ind = %s[i];', cind_name);
                    f.pl('%s[row_ind] = %s[row_ind] + %s[i]*%s[col_ind];', y_name, y_name, val_name, x_name);
                    f.pl('}');
                    
                    %%
                case 'exhaustive_gen'
                    % Write out every non-zero multiplication
                    A = full(A);
                    
                    for r=1:size(A,1)
                        f.p('%s[%i] += 0.0', y_name, r-1)
                        
                        %             if size(A,2) > 50
                        %               [~,ind,val] = find(A(r,:));
                        %               if ~isempty(ind)
                        %                 f.p('+ %.10g*x[%i]', [val;ind-1]);
                        %               end
                        %             else
                        for c=1:size(A,2)
                            if abs(A(r,c)) > 1e-6
                                if A(r,c) == 1
                                    f.p(' + %s[%i]', x_name, c-1)
                                elseif A(r,c) == -1
                                    f.p(' - %s[%i]', x_name, c-1)
                                elseif A(r,c) < 0
                                    f.p(' - %.10g*%s[%i]', abs(A(r,c)), x_name, c-1)
                                elseif A(r,c) > 0
                                    f.p(' + %.10g*%s[%i]', A(r,c), x_name, c-1)
                                end
                            end
                            % end
                        end
                        f.pl(';')
                    end
                    
                case 'FPGA_exhaustive_gen'
                    
                    % Write out every non-zero multiplication
                    A = full(A);
                    
                    for r=1:size(A,1)
                        f.p('%s[%i] += 0.0', y_name, r-1)
                        
                        %             if size(A,2) > 50
                        %               [~,ind,val] = find(A(r,:));
                        %               if ~isempty(ind)
                        %                 f.p('+ %.10g*x[%i]', [val;ind-1]);
                        %               end
                        %             else
                        for c=1:size(A,2)
                            if abs(A(r,c)) > 1e-6
                                if A(r,c) == 1
                                    f.p(' + %s[%i]', x_name, c-1)
                                elseif A(r,c) == -1
                                    f.p(' - %s[%i]', x_name, c-1)
                                elseif A(r,c) < 0
                                    f.p(' - %.10g*%s[%i]', abs(A(r,c)), x_name, c-1)
                                elseif A(r,c) > 0
                                    f.p(' + %.10g*%s[%i]', A(r,c), x_name, c-1)
                                end
                            end
                            % end
                        end
                        f.pl(';')
                    end
                    
                case 'FPGA_matvec_sparse'
                    
                    split_MV_sparse(A, func_name, loop_reset);
                    f.pl('%s_sparse_mv_mult(%s, %s);', func_name, y_name, x_name);
                    
                case 'FPGA_diag'
                    
                    % extract diagonal elements
                    A_diag = diag(A);
                    A_len
                    % store the data
                    val_name = sprintf('%s_diag',Astr);
                    f.add_var(val_name, A_diag, 'type', 'real');
                    
                    % multiply with loop unroll
                    error('User ''FPGA_matvec'' as this feature is yet to be appeared. Not any performancce loss for a rectangular matrix' )
                    
                case 'FPGA_matvec_dense'
                    
                    error('FPGA dense Mat-Vec is yet to be implemented')
            end
        end
    end
end
