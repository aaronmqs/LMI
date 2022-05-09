% Initial configurations
clear; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

% Process Model
Ts = 0.1;

H = tf(1,[1 0.1 2]);

[nH,dH] = tfdata(H,'v');

[A,B,C,D] = tf2ss(nH,dH);

% Condições de desempenho
ta = 5; % settling time
UP = 0.1; % overshoot
xi = sqrt((log(UP))^2/(pi^2 + (log(UP))^2)); % damping coefficient
wn = 4/(xi*ta); % natural frequency
sigma = 4/ta;

% Vertical strip ---beta---sigma---alfa---
alfa = -5*sigma;
beta = -6*sigma;

% LMI (Controller Synthesis)
% A*P + P*A' + B*Z + Z'*B' - 2*h*P < 0
% P > 0
% K = Z/P

% Decision Variables
Z = sdpvar(size(B,2),size(A,1));
P = sdpvar(size(A,2),size(A,2));

Falfa = A*P + P*A' + B*Z + Z'*B' - 2*alfa*P; % Vertical restriction alfa
Fbeta = -A*P - P*A' - B*Z - Z'*B' + 2*beta*P; % Vertical restriction beta

LMI = [ P>=0 ; 
    Falfa <= 0;
    Fbeta <= 0];

optimize(LMI);
checkset(LMI);
Zo = value(Z);
Po = value(P);

% Controller Gain
K = Zo/Po;

M = ss((A+B*K),B,C+D*K,D); % output
M2 = ss((A+B*K),B,K,0); % control signal

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

figure

Gcl = A+B*K;
h1 = alfa;
h2 = beta;
t = -200:0.1:200;
h_1=h1*ones(size(t));
h_2=h2*ones(size(t));

hold on;plot(h_1,t,'r','LineWidth',2); hold on;
hold on;plot(h_2,t,'r','LineWidth',2); hold on;

poles1 = eig(Gcl);
poles_OL1 = eig(A);
scatter(real(poles1),imag(poles1),'rx','LineWidth',2);
hold on
scatter(real(poles_OL1),imag(poles_OL1),'bx','LineWidth',2);
axis([-15 0 -15 15]);
grid on;
title('Pzmap for system 1')











