%% Process parameters

rmp=6.35*10^-3; %[m] raio do motor do carrinho
Kg=3.71; % Reducao das engrenagens do carrinho
Kt=0.00767; %[N.m/A] constante de torque do motor
Rm=2.6; %[ohm] resistencia de armadura do motor
Jm=3.9*10^-7; %[Kg.m^2] momento de inercia do motor
Km=0.00767; %[V.s/rad] forca contra-eletromotriz
Beq=3; %[N.s/m] coeficiente de viscosidade
Kf1=500; %[N/m] rigidez mecanica do primeiro andar
Kf2=500; %[N/m] rigidez mecanica do segundo andar
Mf1=1.160; %[Kg] massa do primeiro andar
Mf2=1.380; %[Kg] massa do segundo andar
Mc=0.65; %[Kg] massa do carrinho

%  dx/dt= A*x + Bu*u + Bw*w
% z = Cz*x + Dzu*u + Dzw*w
% y = Cy*x + Dyu*u + Dyw*w 

[A,Bu,Bw,Cz,Dzu,Dzw,Cy,Dyu,Dyw] = model_2fb(rmp,Kg,Kt,Rm,Jm,Km,Beq,Kf1,Kf2,Mf1,Mf2,Mc);

%% LMI (Controller Synthesis)

% Decision variables
J = sdpvar(size(Bw,2),size(Bw,2));
N = sdpvar(size(Bw,1),size(Bw,1));
M = sdpvar(size(Dzu,2),size(Bw,1));

LMI_H2 = [J >= 0;
          N >= 0;
          [N*A' + M'*Bu' + A*N + Bu*M, N*Cz' + M'*Dzu';
           Cz*N + Dzu*M, -eye(size(Cz,1))] <= 0];

optimize(LMI_H2,trace(J));
checkset(LMI_H2);

J = value(J);
N = value(N);
M = value(M);

K = M/N;

%% H2 plots

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

subplot(2,2,2)
semilogx(wout,Mag1,wout,Mag1_cl,'LineWidth',2);grid on;
axis([10^0 10^3 -150 0])
ylabel("Magnitude(dB)")
str1 = "Open Loop: $||T_{zw}(j \cdot \omega)||_{2}$ = ";
str2 = "Closed Loop: $||T_{zw}(j \cdot \omega)||_{2}$ = ";
title({str1 + "$" + H2 + "$ , " },{str2 + H2_cl},'Interpreter','latex')

subplot(2,2,4)
semilogx(wout,Mag2,wout,Mag2_cl,'LineWidth',2);grid on;
axis([10^0 10^3 -150 0])
ylabel("Magnitude(dB)")

% legend("$||T_{zw}(j \cdot \omega)||$ - Open Loop","$||T_{zw}(j \cdot \omega)||$ - Closed Loop",'FontSize',10,'Interpreter','latex','Location','best')
% sgtitle("$H_{2}$ Control effect on $||T_{zw}(j \cdot \omega)||$",'Interpreter','latex')
xlabel("Frequency $\omega$ (rad/s)",'Interpreter','latex')




















