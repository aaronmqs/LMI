clear, clc, close all;

H = tf(1,[1 -0.1 2]);
[nH,dH] = tfdata(H,'v');

[A,B,C,D] = tf2ss(nH,dH);

Q = 1*eye(2);
R = 1;

% LMI
P = sdpvar(2,2);

F = [ P>=0 ; 
    [(A'*P)+(P*A)+Q P*B; B'*P R]>=0];

optimize(F);
Po = value(P);

K = -inv(R)*B'*value(P);

M = ss((A+B*K),B,C,D);
M2 = ss((A+B*K),B,K,D);

figure
step(M,'r')
hold on
step(M2,'y')