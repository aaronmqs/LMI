% Initial configurations
clear; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
%%

% Process Model
H = tf(1,[1 0.1 2]); % process
[nH,dH] = tfdata(H,'v');

[A,B,C,D] = tf2ss(nH,dH);

% Augmented Matrices
Aa = [A,zeros(size(A,1),size(C,1));-C,zeros(size(C,1),size(C,1))];
Bu = [B;zeros(size(A,1)+size(C,1)-size(B,1),size(B,2))];

%%
% LMI (Controller Synthesis)
% Acl = Aa + Bu*Ka
% Q*Aa' + N'*Bu' + Aa*Q + Bu*N < 0
% Q > 0
% Ka = - N*inv(Q)

% LMI

% Decision Variables
Q = sdpvar(size(Aa,2),size(Aa',1));
N = sdpvar(size(Bu,2),size(Aa,2));

LMI = [ Q>=0 ; 
    Q*Aa' + N'*Bu' + Aa*Q + Bu*N <= 0];

optimize(LMI);
checkset(LMI);
Qo = value(Q);
No = value(N);

% Aumented Gain
Ka = No/Qo;

% Feedback Gain
K = Ka(1:size(A,1));
% Integral gain
H = Ka(size(A,1)+1:end);

%%
% Simulations and plots
Tsim = 30;
simulation = sim('stt_feedback_integral_action');
y = simulation.y;
u = simulation.u;
t = simulation.t;

subplot(2,1,1)
plot(t,y,'LineWidth',3);
legend('y(t)','FontSize',15)
ylabel('Output Signal','FontSize',15)
grid on

subplot(2,1,2)
plot(t,u,'LineWidth',3)
xlabel('Time','FontSize',15)
legend('u(t)','FontSize',15)
ylabel('Control Signal','FontSize',15)
grid on

sgtitle('State Feedback Control - Integral Action')












