classdef splitData < handle
    %
    % Class to write data and functions to files for later use in c-code
    % Harsh changes line 297 to 303 and near line 164
    properties
        % List of strings to be written to the header / c-files
        hfile = {};
        cfile = {};
        
        % List of functions
        functions = {};
        
        % List of data
        data = struct([]);
        
        %         realType = 'double'; % Change to double for double precision
        %         realType = 'float'; % Change to double for double precision
        realType = 'real'; % Change to double for double precision
        fname    = '';
        
        id = 1; % Used to create unique variable names
    end
    
    methods
        function obj = splitData()
        end
        
        function add_var(obj, name, x, varargin)
            % Add variable to be written to the data file
            %
            % add_var(obj, name, x, param, value)
            %
            % name [string]   : name of variable
            % x [double/int]  : data vector / matrix
            %
            % Parameters:
            % type [real]     : {'string', 'int', 'real'}
            % storage_method  : 'dense', 'sparse', 'sparsesuite'
            %                   If not specified, then dense matrices are dense,
            %                   sparse are sparse
            % desc []         : Description of this variable
            %
            % - matrices are stored in column-major form
            
            persistent p
            if isempty(p)
                p = inputParser;
                addRequired(p, 'name',  @ischar);
                addRequired(p, 'x',     @isnumeric);
                expected_types = {'string','int','real'};
                addParameter(p, 'type', 'real', ...
                    @(x) any(validatestring(x,expected_types)));
                expected_storage_methods = {'dense', 'sparse', 'sparsesuite', 'auto'};
                addParameter(p, 'storage_method', 'auto', ...
                    @(x) any(validatestring(x,expected_storage_methods)));
                addParameter(p, 'desc', '', @ischar);
            end
            
            parse(p, name, x, varargin{:});
            [name, x, type, storage_method, desc] = deal(p.Results.name, p.Results.x, p.Results.type, p.Results.storage_method, p.Results.desc);
            
            % Guess the storage method
            if strcmp(storage_method, 'auto')
                if issparse(x)
                    storage_method = 'sparse';
                else
                    storage_method = 'dense';
                end
            end
            
            obj.hl
            
            % Define the sizes
            if ~isempty(desc)
                obj.hl('// %s : %s', name, desc)
            end
            %       obj.hl('const int %s_size[3];', name);
            %       obj.cl('const int %s_size[3] = {%i,%i,%i};', name, size(x,1), size(x,2), numel(x));
            obj.define(sprintf('%s_rows',name), size(x,1), 'int');
            obj.define(sprintf('%s_cols',name), size(x,2), 'int');
            obj.define(sprintf('%s_len', name), length(x(:)), 'int');
            
            switch storage_method
                case 'dense'
                    % Write out a constant dense matrix / vector columnwise
                    i = length(obj.data)+1;
                    obj.data(i).name = name;
                    obj.data(i).x    = x(:);
                    obj.data(i).type = type;
                    
                    obj.alloc(name, type, length(x(:)));
                    
                case 'sparse'
                    % Write out a constant sparse matrix rowwise
                    vec = []; I   = []; len = [];
                    for i = 1:size(x,1)
                        J = find(x(i,:));
                        a = full(x(i,J));
                        
                        len = [len length(J)];
                        I   = [I J];
                        vec = [vec a];
                    end
                    
                    obj.add_var(sprintf('%s_nz_per_row',name), len, 'type', 'int');
                    obj.add_var(sprintf('%s_ind',name), I-1, 'type', 'int');
                    obj.add_var(sprintf('%s_dat',name), vec, 'type', type);
                    
                case 'sparsesuite'
                    error('Have not implemented sparse suite storage method yet')
                    
                otherwise
                    error('Unknown storage method %s', storage_method)
            end
            
        end
        
        function alloc(obj, name, type, len)
            % Allocate memory for a variable
            %
            % name [string]  : name of variable
            % type           : {'string', 'int', 'real'}
            % len            : length of vector
            
            obj.hl('extern %s %s[%i];', obj.type_to_string(type), name, len);
            obj.cl('%s %s[%i];', obj.type_to_string(type), name, len);
        end
        
        function define(obj, name, def, type)
            % Define an item in the header file
            %
            % define(obj, name, def, type)
            %
            %  #define name def
            %
            % type = {'string', 'int', 'real'}
            
            switch type
                case 'char'
                    obj.hl('#define %s "%s"', name, def);
                case 'int'
                    obj.hl('#define %s %i', name, def);
                case 'real'
                    obj.hl('#define %s %f', name, def);
                otherwise
                    error('Unknown type')
            end
        end
        
        function write_to_file(obj, fname)
            obj.fname = fname;
            
            %% Add the data required by all the functions
            for i = 1:length(obj.functions)
                vars = obj.functions{i}.data;
                for j = 1:length(vars)
                    obj.add_var(vars(j).name, vars(j).x, vars(j).args{:});
                end
            end
            
            obj.loadData; % Write a function to load the data
            
            %% Write the binary data file
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%% Temporary solution for FPGA
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            FPGA_head = 1;
            
            if (FPGA_head == 1)
                f = fopen(sprintf('%s.dat', fname), 'w+');
                % struct Header {
                %   unsigned int numVectors;
                %   char info[255];          // Text string describing the file
                % };
                hdr = sprintf(['Optimization data created by SPLIT\n'...
                    'For more information on the SPLIT toolbox visit la.epfl.ch']);
                hdr = {length(obj.data), splitData.format_c_string(hdr, 255)};
                obj.struct_write(f, 'Ic', hdr);
                
                % Write each variable
                for i = 1:length(obj.data)
                    % struct Vector {
                    %   unsigned int type; // 1 => int, 0 => double, 2 => float
                    %   unsigned int len;  // Number of elements
                    %   char name[100];     // Name of the variable
                    %   void *data;      // Vector of sizeof(type)*len
                    % };
                    % typedef struct Vector Vector;
                    switch obj.data(i).type
                        case 'int'
                            dType = 0;
                        case 'real'
                            switch obj.realType
                                case 'float'
                                    dType = 1;
                                case 'double'
                                    dType = 2;
                                case 'real'
                                    dType = 3;
                                otherwise
                                    error('Unkown realType %s', dat.realType)
                            end
                        otherwise
                            error('Unknown type %s', obj.data(i).type)
                    end
                    d = {dType, ...
                        length(obj.data(i).x),...
                        splitData.format_c_string(obj.data(i).name,100),...
                        full(obj.data(i).x)};
                    switch dType
                        case 0
                            obj.struct_write(f, 'IIci', d);
                        case 1
                            obj.struct_write(f, 'IIcf', d);
                        case 2
                            obj.struct_write(f, 'IIcd', d);
                        case 3
                            obj.struct_write(f, 'IIcr', d);
                        otherwise
                            error('Unknown data type')
                    end
                end
                
                fclose(f);
            end
            
            %% Write the header file
            
            f = fopen(sprintf('%s.h', fname), 'w+');
            
            fprintf(f, '/***********************************************\n');
            fprintf(f, ' * Header file automatically generated by SPLIT\n');
            fprintf(f, ' *\n');
            fprintf(f, ' * More information : la.epfl.ch\n');
            fprintf(f, ' ***********************************************/\n');
            
            fprintf(f, '#ifndef __%s_h__\n', fname);
            fprintf(f, '#define __%s_h__\n\n', fname);
            
            fprintf(f, '#include <stdio.h>\n');
            fprintf(f, '#include <stdlib.h>\n');
            fprintf(f, '#include <string.h>\n');
            %   fprintf(f, '#ifdef __APPLE__ \n');
            %  fprintf(f, '#include <Accelerate/Accelerate.h>\n');
            %  fprintf(f, '#endif \n');
            fprintf(f, '#include "user_matrix_ops.h"\n');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%% Temporary solution for FPGA
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %fprintf(f, '#include "user_splitLoad.h"\n');
            %fprintf(f, '#include "user_probdata_FPGA.h"\n');
            
            %  fprintf(f, '#include "user_ldl.h"\n\n');
            
            
            
            
            fprintf(f, '#define real %s\n', obj.realType);
            fprintf(f, '#define REAL %s\n', obj.realType);
            
            
            %fprintf(f, 'typedef %s REAL;\n', obj.realType);
            switch obj.realType
                case 'double'
                    fprintf(f, '#define printVec printDoubleVec\n');
                case 'float'
                    fprintf(f, '#define printVec printFloatVec\n');
            end
            
            
            %       fprintf(f, '#define SIZE(var, dim) var##_size[dim-1]\n');
            %       fprintf(f, '#define NUMEL(var) var##_size[2]\n');
            
            for i = 1:length(obj.hfile)
                fprintf(f, obj.hfile{i}{:});
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%% Temporary solution for declaring the variable
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fprintf(f, '#ifdef lapack_linsolve \n');
            fprintf(f, 'extern __CLPK_integer ipiv[nn_lp];\n');
            fprintf(f, '#endif \n');
            fprintf(f, 'extern double *Lx_ss;\n');
            fprintf(f, 'extern int *Li_ss;\n');
            
            % Print all the function prototypes
            fprintf(f, '\n\n// Functions\n');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%% Temporary solution for FPGA
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            FPGA_head = 1;
            if (FPGA_head == 1)
                
                for i = 1:length(obj.functions)-1
                    fprintf(f, obj.functions{i}.print_hfile);
                end
            else
                
                for i = 1:length(obj.functions)
                    fprintf(f, obj.functions{i}.print_hfile);
                end
                
            end
            fprintf(f, '\n\n');
            
            fprintf(f, '#endif\n');
            
            fclose(f);
            
            %% Write the c-file
            
            f = fopen(sprintf('%s.c', fname), 'w+');
            
            fprintf(f, '/***********************************************\n');
            fprintf(f, ' * Header file automatically generated by SPLIT\n');
            fprintf(f, ' *\n');
            fprintf(f, ' * More information : la.epfl.ch\n');
            fprintf(f, ' ***********************************************/\n');
            
            fprintf(f, '#include "%s.h"\n\n', fname);
            
            
            %%%% We do not need the following if we write data in c file.
            
%             for i = 1:length(obj.cfile)
%                 fprintf(f, obj.cfile{i}{:});
%             end
            
            
            fprintf(f, '\n\n');
            
            % Print all the functions
            FPGA_head = 1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%% Temporary solution for FPGA not to print load data - old
            %%%%%% split systems
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (FPGA_head == 1)
                for i = 1:length(obj.functions)-1
                    fprintf(f, obj.functions{i}.print_cfile);
                end
                
            else
                for i = 1:length(obj.functions)
                    fprintf(f, obj.functions{i}.print_cfile);
                end
            end
            
            write_FPGA(obj,f)
            
            fclose(f);
            
        end
        
        function add_function(dat, f)
            for i = 1:length(f)
                funcs = f(i).getFunctions;
                dat.functions = {dat.functions{:} funcs{:}};
            end
        end
    end
    
    methods (Static, Hidden = true, Access = 'private')
        function str = format_c_string(str, len)
            i = length(str);
            str = sprintf('%-*s',len, str);
            if i < len
                str(i+1) = 0;
            else
                str(len+1) = 0;
            end
        end
    end
    
    methods (Access = 'private', Hidden = true)
        function name = uniqueVarName(dat)
            name = sprintf('_var_%i_', dat.id);
            dat.id = dat.id + 1;
        end
        
        function loadData(obj)
            % Write out a function load all data into the memory
            
            f = coderFunc('void loadData()');
            f.pl('Data *dat = splitLoad("%s.dat");', obj.fname);
            for i = 1:length(obj.data)
                f.pl('copyVar(%s, dat, "%s");', obj.data(i).name, obj.data(i).name);
            end
            f.pl('freeData(dat);')
            
            obj.add_function(f);
        end
        
        function str = type_to_string(obj, type)
            switch type
                case 'string'
                    str = 'char';
                case 'int'
                    str = 'int';
                case 'real'
                    str = obj.realType;
                otherwise
                    error('Unknown type %s', type)
            end
        end
    end
    
    methods
        function h(obj, varargin)
            % Print to the header file
            if isempty(varargin), return; end
            obj.hfile{end+1} = varargin;
        end
        function hl(obj, varargin)
            % Print a line to the header file
            if length(varargin) > 0
                obj.h(varargin{:})
            end
            obj.h('\n')
        end
        function c(obj, varargin)
            % Print to the c-file
            if isempty(varargin), return; end
            obj.cfile{end+1} = varargin;
        end
        function cl(obj, varargin)
            % Print a line to the c-file
            if length(varargin) > 0
                obj.c(varargin{:})
            end
            obj.c('\n')
        end
    end
end


