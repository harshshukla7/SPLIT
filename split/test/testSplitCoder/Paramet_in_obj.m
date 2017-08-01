%% Lets say we want to solve the following optimization problem:
%%%%%%%%%% minimize (z'z + r'z) 
%%%%%%%%%% subject to Az <= b 
%%%%%%%%%% where optimization variable is z
%%%%%%%%%% and parameters are in r

clear all 
close all
n = 4; %number of optimization variables
A = rand(n,n);
b = rand(n,1);

switch cas
    
    case 1
r = parameter(2*n,1); % this parameter will appear affine in the objective function.
x = splitvar(n,1);
u = splitvar(n,1);

objQ1 =  (x .* x);
objQ2 =  (u .* u);
objA =  2*(r(1:n) .* x); 
obj = sum(objQ1(:)) + sum(objQ2(:)) + sum(objA(:));

A*x <= b;
        
minimize(obj);
        
% Generate the problem data
prob = splitProb.genProblem;

    case 2