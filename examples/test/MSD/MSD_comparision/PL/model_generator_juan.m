function [model_c] = model_generator_juan(N_carts)
% this function accepts the number of carts (masses) and
% returns continious state space model of a mass spring damper system

design.n_states = N_carts*2;
design.m_inputs = N_carts;

k = 1*ones(1,N_carts); % spring constants
b = 0*0.2*ones(1,N_carts); % damping constants
m = 1*ones(1,N_carts); % damping constants


k(N_carts+1) = 0; % the last cart is connected only to previous one
b(N_carts+1) = 0;

% state matrix generation
A = kron(diag(ones(1,N_carts)),[0 1; 0 0]); 

for i=1:N_carts
    A(i*2,i*2-1) = -(k(i)+k(i+1))/m(i);
    A(i*2,i*2) = -(b(i)+b(i+1))/m(i); 
    
    if(i>1)
        A(i*2,i*2-3) = k(i)/m(i);
        A(i*2,i*2-2) = b(i)/m(i);
    end
    
    if(i<N_carts)
        A(i*2,i*2+1) = k(i+1)/m(i);
        A(i*2,i*2+2) = b(i+1)/m(i);
    end
    
    
end

% input matrix generation
B = zeros(N_carts*2,N_carts);
for i=1:N_carts
    B(i*2,i) = 1 / m(i);
end

% output matrix generation
C=eye(N_carts*2);

% Direct transition matrix generation
D = zeros(N_carts*2, N_carts);

model_c = ss(A,B,C,D);

end
