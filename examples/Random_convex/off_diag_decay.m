
clear all
close all

band = 300;
alpha = 1e-2;
steps = 30;
th = 1e-6;

n =  60;
m =  30;
N =  70;


[A,B,C,D] = ssdata(drss(n,n,m));

x = splitvar(n,N);
u = splitvar(m,N-1);

-10 <= x(:) <= 10;
-1 <= u(:) <= 1;

x0 = parameter(n,1);

%r = parameter(n,N);
x(:,1) == x0;
x(:,2:end) == A*x(:,1:end-1) + B*u;

Q = eye(n);
R = eye(m);
% objX = x .* (Q*x);
% objU = u .* (R*u);
% obj = sum(objX(:)) + sum(objU(:));

%         objX = x.*(Q*(x));
%         objU = u.*((R*u));
%         obj = sum(objX(:)) + sum(objU(:));
obj = 0;
for i=1:N-1
    obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i);
end
obj = obj + x(:,N)'*Q*x(:,N);


minimize(obj);

prob = splitProb.genProblem;

dat = prob.coderData;
K1 = sparse(dat.Q+dat.L'*dat.L);
K = sparse([K1 dat.A'; dat.A zeros(size(dat.A,1))]);

kp = symrcm(K);
KK = K(kp,kp);

tic
[LL,DD,PP] = ldl(KK);
tp3_time = toc;


b = rand(size(KK,1),1);

x1 = PP'\(LL'\(DD\(LL\(PP\(b)))));




tic
for i=0:1000
x1 = PP'\(LL'\(DD\(LL\(PP\(b)))));
end
tp1_time = toc;

%% 


% KK_in = speye(size(KK,1))/KK;
% KK_in_band = KK_in.*Band_mat;
% x3 = KK_in_band*b;
% norm(x3-x1)
% x2 = KK_in*b; norm(x1-x2);

% tic
% for i=0:1000
% x3 = KK_in_band*b;
% end
% tp2_time = toc;


%% compute the inverse using Newton iteration

Band_mat = createbandOPT(size(K,1),band);
IK = 2*speye(size(K,1));
KK_in_newton = alpha*KK;
KK_norm = normest(KK);

err=zeros(steps,1);

tic
for i=1 : steps 
     
    KK_in_newton = Band_mat.*(KK_in_newton*(IK-KK*KK_in_newton));
    KK_in_newton(KK_in_newton<th)=0;
    
    %err(i) = normest(IK-2*KKKK)/KK_norm;
        
end
tp4_time = toc;


KK_in_newton =  sparse(KK_in_newton);
x4=KK_in_newton*b;
relative_norm = norm(x4-x1)/norm(x1)

 tic
 for i=0:1000
 x4 = KK_in_newton*b;
 end
 tp2_time = toc;

%%
% figure
% Kinv = eye(size(K,1))/KK;
% surf((Kinv))
% title('AL-KKT')


% K = [dat.Q dat.A'; dat.A zeros(size(dat.A,1))];
% kp = symrcm(K);
% KK = K(kp,kp);

% figure
% Kinv = eye(size(K,1))/KK;
% surf((Kinv))
% title('L-KKT');