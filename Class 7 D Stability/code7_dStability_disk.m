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

% Disk of radius r
r = 0.2*sigma;
q = sigma;

% LMI (Controller Synthesis)
% A*P + P*A' + B*Z + Z'*B' - 2*h*P < 0
% P > 0
% K = Z/P

% Decision Variables
P = sdpvar(size(A,2),size(A,2));
Z = sdpvar(size(B,2),size(A,2));

LMI = [ P>=0 ; 
    [-r*P,q*P+A*P+B*Z;
    q*P+P*A'+Z'*B',-r*P] <= 0];

optimize(LMI);
checkset(LMI);
Zo = value(Z);
Po = value(P);

% Controller Gain
K = Zo/Po;

M = ss((A+B*K),B,C+D*K,D,Ts); % output
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

figure
pzmap(H)
hold on
pzmap(M)
hold on
circle(-q,0,r)

function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold off
end




