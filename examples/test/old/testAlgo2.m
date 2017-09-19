clear all

algo = splitAlgo;

n = 5;

x = algoVar(n, 1, 'x');
y = algoVar(n, 1, 'y');
z = algoVar(1, 1, 'z');

M = algoConst(randn(n,2*n),'M');
c = algoConst(randn(n,1),'c');
Q = randn(n); 
Q = algoConst(Q*Q','Q');

algo.assign(x, norm(M*[x;y] + 3));

algo.comment('What an amazing comment!');
algo.assign(y, ones(n,1)*norm(x,inf));
algo.assign(z, x'*Q*x + y'*c);

i = algoVar(1,1,'i');
algo.whileLoop(i >= 0);
  algo.comment('Here''s a comment!');
  algo.assign(x, 3+i);
algo.endIf;

op = algo.ops{1};
algo.flatten(op)
