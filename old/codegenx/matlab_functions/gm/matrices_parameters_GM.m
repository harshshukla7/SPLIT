function GM = matrices_paramaters_GM(dat)
%% build matrices for fast gradient method
% matrices to build - reference: see Stefan Richter dissertation thesis,
% chapter 8 - Certification of input-constrained MPC
%             Hg  = Bg'Qg'Bg + Rg
%             M1g = Ag'Qg Bg
%             M2g = Ag'Qg Ag

A_g = [];
B_g = zeros((dat.N)*size(dat.B,1),(dat.N-1)*size(dat.B,2));
Q_g = kron(dat.Q,eye(dat.N));
R_g = kron(dat.R,eye(dat.N-1));

% temporary matrix needed to construct the matrix B_g
B_help = [];
for j = 1:dat.N-1
        B_help = [dat.A^(j-1)*dat.B, B_help];
end

% construction of matrices Ag and Bg
for i = 1: dat.N 
    A_g = [A_g' dat.A^(i-1)']';
    
    if i > 1
        B_g(1+(i-1)*size(dat.B,1):(i)*size(dat.B,1),1:(i-1)*size(dat.B,2)) = B_help(:,end-(i-1)*size(dat.B,2)+1:end);
    end
end

% final matrices
GM.H_g = B_g'*Q_g*B_g+R_g;
GM.M1_g = A_g'*Q_g*B_g;
GM.M2_g = A_g'*Q_g*A_g;
GM.mu = min(eig(GM.H_g)); % convexity paramater
GM.L = max(eig(GM.H_g));  % lipschitz constant
GM.alpha0 = (1-sqrt(GM.mu/GM.L))/5*4; 
GM.epsilon = dat.epsilon;
GM.step = dat.step;
GM.horizon = dat.N-1;
GM.num_inputs = dat.m;
GM.init_conditions = dat.x0;

end
