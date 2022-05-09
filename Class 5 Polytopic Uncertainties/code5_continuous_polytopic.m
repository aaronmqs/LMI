clear, clc, close all;

delta = -2*rand + 1; % random number in (-1,1)
A1 = [-4 4;-5 0] + [-2 2;-1 4];
A2 = [-4 4;-5 0] - [-2 2;-1 4];
A3 = [-4 4;-5 0] + delta*[-2 2;-1 4]; % test matrix
B = [1;2];
C = [1 1];
D = 0;

% LMI1 - Analysis
P = sdpvar(2,2);

LMI0 = P>=0 ; 
LMI1 = A1*P+P*A1<=0;
LMI2 = A2*P+P*A2<=0;

LMI_analysis = [LMI0,LMI1,LMI2];

% LMI_analysis = [ P>=0 ; 
%     A1*P+P*A1<=0;
%     A2*P+P*A2<=0];

optimize(LMI_analysis);

checkset(LMI_analysis);

% LMI1 - Synthesis
Q = sdpvar(2,2);
Y = sdpvar(1,2);

LMI0 = Q>=0 ;

LMI1 = Q*A1'+A1*Q+Y'*B'+B*Y <= 0;

LMI2 = Q*A2'+A2*Q+Y'*B'+B*Y <= 0;

LMI_synthesis = [LMI0,LMI1,LMI2];

% LMI_synthesis = [ Q>=0 ; 
%     Q*A1'+A1*Q+Y'*B'+B*Y <= 0;
%     Q*A2'+A2*Q+Y'*B'+B*Y <= 0;];

optimize(LMI_synthesis)

Qo = value(Q);
Yo = value(Y);

K = - Yo/Qo;

M1 = ss((A1-B*K),B,C,D); % output
Mu1 = ss((A1-B*K),B,K,D); % control signal

M2 = ss((A2-B*K),B,C,D); % output
Mu2 = ss((A2-B*K),B,K,D); % control signal

M3 = ss((A3-B*K),B,C,D); % output
Mu3 = ss((A3-B*K),B,K,D); % control signal

subplot(2,1,1)
step(M1,'r')
title('Closed Loop Output - A1')
set(findall(gcf,'type','line'),'linewidth',3);
grid on

subplot(2,1,2)
step(Mu1,'y')
title('Control Signal')
set(findall(gcf,'type','line'),'linewidth',3);
grid on

figure

subplot(2,1,1)
step(M2,'r')
title('Closed Loop Output - A2')
set(findall(gcf,'type','line'),'linewidth',3);
grid on

subplot(2,1,2)
step(Mu2,'y')
title('Control Signal')
set(findall(gcf,'type','line'),'linewidth',3);
grid on

figure

subplot(2,1,1)
step(M3,'r')
title('Closed Loop Output - Random Matrix between A1-A2')
set(findall(gcf,'type','line'),'linewidth',3);
grid on

subplot(2,1,2)
step(Mu3,'y')
title('Control Signal')
set(findall(gcf,'type','line'),'linewidth',3);
grid on


