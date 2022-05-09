clear, clc, close all;

Ts = 0.1;

H = tf(1,[1 -1 2]);
Hd = c2d(H,Ts);

[nH,dH] = tfdata(Hd,'v');

[A,B,C,D] = tf2ss(nH,dH);

% LMI (síntese de controlador)
% [-Q,Q*A'+X*B';A*Q+B*X,-Q] <= 0
% Q >= 0
% K = X/Q

% variáveis de decisão
Q = sdpvar(2,2);
N = sdpvar(1,2);

LMI = [ Q>=0 ; 
    [-Q Q*A'+N'*B';A*Q+B*N -Q] <= 0];

optimize(LMI)
Qo = value(Q);
No = value(N);

K = No/Qo;

checkset(LMI)

M = ss((A+B*K),B,C,D,Ts); % output
M2 = ss((A+B*K),B,K,0,Ts); % control signal

figure
step(M,'r')
title('Output')
set(findall(gcf,'type','line'),'linewidth',3);
grid on

figure
step(M2,'y')
title('Control Signal')
set(findall(gcf,'type','line'),'linewidth',3);
grid on

figure
step(H)
title('Open Loop Output')
set(findall(gcf,'type','line'),'linewidth',3);
grid on

