%% Riccati recursion methods:

function [K_cont,P] = riccati_recurison_quad(A,B,P,Q,R,S,N)
        
            K_cont = cell(1,N);
       
        for i=N-1:-1:0
            R = R + B'*P*B;
            K = -(R\(S+B'*P*A));
            P = Q + A'*P*A - K'*R*K;
            K_cont{1,i+1} = K;
        end   
        
        
        
%         for i=0:N-1
%             
%             u = K_cont{1,i+1}*x;
%             x = A*x + B*u;
%     
%           
%         end
end
    
    
