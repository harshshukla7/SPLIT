clear 

% Create an object to store data that will be written to the file
% fname.dat in binary format and static memory allocations to fname.h
dat = splitData('splitData1');

% When you call dat.add, it will write a line to the file fHdr
%   double x[10];
% or
%   int x[10];

% Add a bunch of *vectors*, giving them names
dat.add('x',randn(10,1));
dat.add('sol',randn(20,1));

% Write matrix variables
A = randn(10,3);
dat.add('Arow', vec(A')); % Row major format
dat.add('Acol', vec(A));  % Column major format

% Write integer vectors
I = ceil(100*rand(30,1));
dat.add('I', I, 'int');

% Write the data to file
dat.writeFile;
