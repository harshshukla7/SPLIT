%% Riccati recursion methods:

function [P] = riccati_recurison_exended(A,B,P,Q,R,S,N,b,p,q,s)


    K_cont = cell(1,N);
    k_cont = cell(1,N);   
    
for i=N-1:-1:0
    
    R = R + B'*P*B;
    K = -(R\(S + B'*P*A));
    P = Q + A'*P*A - K'*R*K;
    k = -(R\(s + B'*(P*b + p)));
    p = q + A'*(P*b + p) - K'*R*k;
    
    K_cont{1,i+1} = K;
    k_cont{1,i+1} = k;
    
end


for i=0:N-1
    
    u = K_cont{1,i+1}*x + k_cont{1,i+1};
    x = A*x + B*u + b;
    
end


