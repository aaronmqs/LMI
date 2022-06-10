clc; clear; close all;

% Compares the separeted effects of H2 and Hinf Control
TwoFloorBuilding_Hinf;
TwoFloorBuilding_H2;
legend("$||T_{zw}(j \cdot \omega)||$ - Open Loop","$||T_{zw}(j \cdot \omega)||$ - Closed Loop",'FontSize',10,'Interpreter','latex','Location','best')

%% H2Hinf Control

% Decision variables
Q = sdpvar(size(A,1),size(A,2));
N = sdpvar(size(Bu,2),size(A,2));
J = sdpvar(size(Bw,2),size(Bw,2));
delta = sdpvar(1,1);

% LMI
LMI_Hinf = [Q >= 0;
            [Q*A' + N'*Bu' + A*Q + Bu*N, Bw, Q*Cz' + N'*Dzu';
             Bw', -delta*eye(size(Bw,2),size(Bw,2)), Dzw';
             Cz*Q + Dzu*N, Dzw, -eye(size(Dzw,1),size(Dzw,1))] <= 0];

LMI_H2 = [J >= 0;
          Q >= 0;
          [Q*A' + N'*Bu' + A*Q + Bu*N, Q*Cz' + N'*Dzu';
           Cz*Q + Dzu*N, -eye(size(Cz,1))] <= 0];

LMI_H2Hinf = [LMI_Hinf;LMI_H2];

optimize(LMI_H2Hinf,trace(J));
checkset(LMI_H2Hinf);

Q = value(Q);
N = value(N);

% Controller gain
K = N/Q;

%% H2Hinf plots

w = 1:0.01:10^3;

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

H2 = sqrt(trace(covar(Tzw,1))); % open loop H2 norm;
if H2 == Inf
    H2 = "\infty";
end
H2_cl = sqrt(trace(covar(Tzw_cl,1))); % closed loop H2 norm;
Hinf = hinfnorm(Tzw); % Open Loop Hinf norm
if Hinf == Inf
    Hinf = "\infty";
end
Hinf_cl = hinfnorm(Tzw_cl); % Closed Loop Hinf norm

figure

subplot(2,1,1)
semilogx(wout,Mag1,wout,Mag1_cl,'LineWidth',2);grid on;
axis([10^0 10^3 -150 0])
ylabel("Magnitude(dB)")
str1 = "Open Loop: $||T_{zw}(j \cdot \omega)||_{\infty}$ = ";
str2 = "Closed Loop: $||T_{zw}(j \cdot \omega)||_{\infty}$ = ";
str3 = "Open Loop: $||T_{zw}(j \cdot \omega)||_{2}$ = ";
str4 = "Closed Loop: $||T_{zw}(j \cdot \omega)||_{2}$ = ";
title({str1 + "$" + Hinf + "$" ...
    str2 + Hinf_cl ...
    str3 + "$" + H2 + "$" ...
    str4 + H2_cl},'Interpreter','latex')

subplot(2,1,2)
semilogx(wout,Mag2,wout,Mag2_cl,'LineWidth',2);grid on;
axis([10^0 10^3 -150 0])
ylabel("Magnitude(dB)")

legend("$||T_{zw}(j \cdot \omega)||$ - Open Loop","$||T_{zw}(j \cdot \omega)||$ - Closed Loop",'FontSize',10,'Interpreter','latex','Location','best')
xlabel("Frequency $\omega$ (rad/s)",'Interpreter','latex')