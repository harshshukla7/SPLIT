clear all

dat = splitData('splitData');
global f
f = dat.fHdr;
gen = splitGen(dat);

A = randn(10); A = A*A';
A = triu(tril(A,3),-3);

gen.gen_solve_banded(A, 'A', 'test_banded');

b = randn(size(A,1),1);
dat.add('b', b);

dat.add('sol', A\b);
dat.define('n',size(A,1),'int');

dat.writeFile;
