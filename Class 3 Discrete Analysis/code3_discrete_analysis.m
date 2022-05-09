clear, clc, close all;

Ts = 0.1;

H = tf(1,[1 -1 2]);
Hd = c2d(H,Ts);

[nH,dH] = tfdata(Hd,'v');

[A,B,C,D] = tf2ss(nH,dH);

% LMI

% Analysis
% variáveis de decisão
P = sdpvar(2,2);

LMI = [ P>=0 ; 
    A'*P*A-P<=0];

optimize(LMI)
Po = value(P)

eig(Po)
checkset(LMI)











