clear

% Create an object to store data that will be written to the file
% fname.dat in binary format and static memory allocations to fname.h
dat = splitData;

%%
% Generate a matrix-vector multiplication
A = sprandn(100,44,0.1);
b = randn(100,1);
x = randn(44,1);
dat.add_var('x', x);
dat.add_var('sol', A*x+b);

% We create functions that are derived from coderFunc, and then add them to
% the splitData object
f(1) = coderFunc_times('func_blas',       A, 'bConst', b, 'method', 'blas');
f(2) = coderFunc_times('func_forloop',    A, 'bConst', b, 'method', 'for_loops');
f(3) = coderFunc_times('func_sparse',     A, 'bConst', b, 'method', 'sparse');
f(4) = coderFunc_times('func_exhaustive', A, 'bConst', b, 'method', 'exhaustive');
dat.add_function(f);

% You can also #define anything this way
dat.define('n', size(A,1), 'int');

%% Test solution of linear systems of equations
A = sprandn(100,100,0.1); A = A + A';

% Store the test vector and solution
x = randn(100,1);
dat.add_var('x_lin', x);
dat.add_var('y_sol', A\x);

clear f
f(1) = coderFunc_mldivide('func_solve_inv', A, 'method', 'invert');
f(2) = coderFunc_mldivide('func_solve_ldl', A, 'method', 'ldl');

dat.add_function(f);

%% Here's a bunch of examples on how to add variables to a data file
% add_var a bunch of vectors & matrices, giving them names
dat.add_var('aaa',randn(10,1),  'type', 'real', 'desc', 'This is a vector aaa');
dat.add_var('bbb',randn(20,1),'type', 'real');

% Write matrix variables
A = randn(10,3);
dat.add_var('Arow', A','type', 'real'); % Row major format
dat.add_var('Acol', A, 'type', 'real'); % Column major format

% Write sparse matrix
A = sprandn(50,50,0.2);
dat.add_var('Asparse', A,'type', 'real', 'desc', 'This is a sparse matrix');  % Column major format

% Write integer vectors
I = ceil(100*rand(30,1));
dat.add_var('I', I, 'type', 'int');
% Force sparse for a dense matrix
dat.add_var('Isparse', I, 'type', 'int', 'storage_method', 'sparse', 'desc', 'This is a dense integer matrix stored in sparse format');


%% Write the data to file
dat.write_to_file('test');

%% Compile
system('gcc -o test -framework Accelerate testSplitData.c test.c splitLoad.c');
system('./test');
