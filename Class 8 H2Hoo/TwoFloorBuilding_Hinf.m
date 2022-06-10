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

%% Hinf plots

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
Hinf = hinfnorm(Tzw); % Open Loop Hinf norm
if Hinf == Inf
    Hinf = "\infty";
end
Hinf_cl = hinfnorm(Tzw_cl); % Closed Loop Hinf norm

subplot(2,2,1)
semilogx(wout,Mag1,wout,Mag1_cl,'LineWidth',2);grid on;
axis([10^0 10^3 -150 0])
ylabel("Magnitude(dB)")
str1 = "Open Loop: $||T_{zw}(j \cdot \omega)||_{\infty}$ = ";
str2 = "Closed Loop: $||T_{zw}(j \cdot \omega)||_{\infty}$ = ";
title({str1 + "$" + Hinf + "$ , " },{str2 + Hinf_cl},'Interpreter','latex')

subplot(2,2,3)
semilogx(wout,Mag2,wout,Mag2_cl,'LineWidth',2);grid on;
axis([10^0 10^3 -150 0])
ylabel("Magnitude(dB)")

% legend("$||T_{zw}(j \cdot \omega)||$ - Open Loop","$||T_{zw}(j \cdot \omega)||$ - Closed Loop",'FontSize',10,'Interpreter','latex','Location','best')
% sgtitle("$H_{\inf}$ Control effect on $||T_{zw}(j \cdot \omega)||$",'Interpreter','latex')
xlabel("Frequency $\omega$ (rad/s)",'Interpreter','latex')





















