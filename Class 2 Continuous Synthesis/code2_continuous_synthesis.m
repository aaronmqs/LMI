clear
clc;
close all;

H = tf(1,[1 0.1 2]); % process
[nH,dH] = tfdata(H,'v');

[A,B,C,D] = tf2ss(nH,dH);

% LMI (síntese de controlador)
% Q*A' + A*Q + Y'B' + B*Y < 0
% Q > 0
% K = - Y*inv(Q)

% LMI

% variáveis de decisão
Q = sdpvar(2,2);
Y = sdpvar(1,2);

LMI = [ Q>=0 ; 
    Q*A'+A*Q+Y'*B'+B*Y<=0];

optimize(LMI);
Qo = value(Q);
Yo = value(Y);

K = - Yo/Qo;

M = ss((A-B*K),B,C,D); % output
%M = ss((A-B*K),B,(C-D*K),D); % output
M2 = ss((A-B*K),B,K,D); % control signal

figure
step(M,'r')
hold on
step(M2,'y')
legend('Output','Control Signal')

figure
step(H)
