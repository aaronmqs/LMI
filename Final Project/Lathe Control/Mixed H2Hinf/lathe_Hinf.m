%% Process parameters
close all
Mf = 0.1; % Kg
Mp = 1; % Kg
Kf = 500; % N/m
Kp = 10; % N/m
Cf = 0.0001; % N.s/m
Cp = 0.1; % N.s/m


A = [0 1 0 0;
    -Kf/Mf -Cf/Mf Kf/Mf Cf/Mf;
    0 0 0 1;
    Kf/Mp Cf/Mp -(Kp+Kf)/Mp -(Cp+Cf)/Mp];

Bu = [0 1/Mf 0 0]';

Bw = [0 1/Mf 0 0]';

Cz = [1 0 0 0;
    0 0 1 0];

Dzu = [0 0]';

Dzw = [0 0]';

%  dx/dt= A*x + Bu*u + Bw*w
% z = Cz*x + Dzu*u + Dzw*w
% y = Cy*x + Dyu*u + Dyw*w 

%% LMI (Controller Synthesis)

% Decision variables
Q = sdpvar(size(A,1),size(A,2));
N = sdpvar(size(Bu,2),size(A,2));
delta = sdpvar(1,1);

% LMI
LMI_Hinf = [Q >= 0;
            [Q*A' + N'*Bu' + A*Q + Bu*N, Bw, Q*Cz' + N'*Dzu';
             Bw', -delta*eye(size(Bw,2),size(Bw,2)), Dzw';
             Cz*Q + Dzu*N, Dzw, -eye(size(Dzw,1),size(Dzw,1))] <= 0];

optimize(LMI_Hinf,delta);
checkset(LMI_Hinf);

Q = value(Q);
N = value(N);

% Controller gain
K = N/Q;

gama = sqrt(double(delta));

%% Hinf plots

w = 1:0.001:10^3;

% Tzw - Open Loop
Tzw = ss(A,Bw,Cz,Dzw);
[mag,~,~] = bode(Tzw,w);
Mag1=squeeze(20*log10(mag(1,1,:)));
Mag2=squeeze(20*log10(mag(2,1,:)));


% Tzw - Closed Loop
Tzw_cl = ss(A+Bu*K,Bw,Cz+Dzu*K,Dzw);
[mag,~,wout] = bode(Tzw_cl,w);
Mag1_cl=squeeze(20*log10(mag(1,1,:)));
Mag2_cl=squeeze(20*log10(mag(2,1,:)));
Hinf = hinfnorm(Tzw); % Open Loop Hinf norm
if Hinf == Inf
    Hinf = "\infty";
end
Hinf_cl = hinfnorm(Tzw_cl); % Closed Loop Hinf norm

subplot(2,2,1)
semilogx(wout,Mag1,wout,Mag1_cl,'LineWidth',2);grid on;
axis([10^0 10^3 -150 50])
ylabel("Magnitude(dB)")
str1 = "Open Loop: $||T_{zw}(j \cdot \omega)||_{\infty}$ = ";
str2 = "Closed Loop: $||T_{zw}(j \cdot \omega)||_{\infty}$ = ";
title({str1 + "$" + Hinf + "$ , " },{str2 + Hinf_cl},'Interpreter','latex')

subplot(2,2,3)
semilogx(wout,Mag2,wout,Mag2_cl,'LineWidth',2);grid on;
axis([10^0 10^3 -150 50])
ylabel("Magnitude(dB)")

% legend("$||T_{zw}(j \cdot \omega)||$ - Open Loop","$||T_{zw}(j \cdot \omega)||$ - Closed Loop",'FontSize',10,'Interpreter','latex','Location','best')
% sgtitle("$H_{\inf}$ Control effect on $||T_{zw}(j \cdot \omega)||$",'Interpreter','latex')
xlabel("Frequency $\omega$ (rad/s)",'Interpreter','latex')





















