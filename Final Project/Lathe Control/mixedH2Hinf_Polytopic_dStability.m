%% Initial configurations
clear; clc; close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex'); % interpretador LateX para os eixos
set(groot, 'defaultLegendInterpreter','latex'); % interpretador LateX para a legenda
set(groot,'defaultTextInterpreter','latex'); % interpretador de texto padrao (LateX)

%% Process Model
% Mf = 0.1; % Kg
% Mp = 1; % Kg
% Kf = 500; % N/m
% Kp = 10; % N/m
% Cf = 0.0001; % N.s/m
% Cp = 0.1; % N.s/m

Mf = 0.1; % Kg
Mp = 1 + (0.2*rand(1,50) - 0.1); % Kg
% Mp = 1;
Kf = 500; % N/m
Kp = 10 + (2*rand(1,50) - 1); % N/m
% Kp = 500;
Cf = 0.0001; % N.s/m
Cp = 0.1; % N.s/m
n = 2; % number of uncertainties

H = cell(n); % Polytope Vertices
Mp_i = H; % uncertainties in the damping coefficient
Kp_i = H; % uncertainties in the static gain
poles_OL = H; % Open Loop Poles
% State Space Matrices
A = H;
B = H;
Bw = H;
C = H;
D = H;

% Polytop plotting
figure
scatter(Mp,Kp,'filled');
ylabel("Spring Constant (N.m)")
xlabel("Piece Mass (Kg)")
title("Polytopic uncertainties")
hold on
grid on

for i = 1:n
    for j = 1:n
        if mod(j,2) == 0
            Mp_i{i,j} = max(Mp);
        else 
            Mp_i{i,j} = min(Mp);
        end
        if mod(i,2) == 0
            Kp_i{i,j} = min(Kp);
        else 
            Kp_i{i,j} = max(Kp);
        end
        scatter(Mp_i{i,j},Kp_i{i,j},'rx','LineWidth',2);
        hold on
        A{i,j} = [0 1 0 0;
                -Kf/Mf -Cf/Mf Kf/Mf Cf/Mf;
                0 0 0 1;
                Kf/Mp_i{i,j} Cf/Mp_i{i,j} -(Kp_i{i,j}+Kf)/Mp_i{i,j} -(Cp+Cf)/Mp_i{i,j}];

        B{i,j} = [0 1/Mf 0 0]';
        Bw{i,j} = [0 1/Mf 0 0]';
        C{i,j} = [1 0 0 0];
        D{i,j} = 0;
        poles_OL{i,j} = eig(A{i,j});
    end
end

%% Performance parameters
t_min = 0.1; % min settling time
t_max = 500; % max settling time
UP = 0.05; % overshoot
zeta = sqrt((log(UP))^2/(pi^2 + (log(UP))^2)); % damping coefficient
h1 = -4/t_max; % Vertical Strip
h2 = -4/t_min; % Vertical Strip
theta = acos(zeta); % Conic Sector angle
tmax = 45; % max value of vector t

%% LMI (Controller Synthesis)

% Decision Variables
P = sdpvar(size(A{1,1},1),size(A{1,1},1));
Z = sdpvar(size(B{1,1},2),size(A{1,1},1));
J = sdpvar(size(Bw{1,1},2),size(Bw{1,1},2));
% delta = sdpvar(1,1);
delta = 1;

LMIJ = J >= 0;

LMI0 = P >= 0;

LMI1 = cell(n); % Matrix of LMIs
LMI_Hinf = LMI1;
LMI_H2 = LMI1;
LMI_H2Hinf = LMI1;


for i = 1:n
    for j = 1:n

        % Conic Sector Constraints
        % Left Vertical Strip Constraint
        % Right Vertical Strip Constraint
        LMI1{i,j} = [[sin(theta)*(A{i,j}*P + B{i,j}*Z + P*A{i,j}' + Z'*B{i,j}'),cos(theta)*(A{i,j}*P + B{i,j}*Z - P*A{i,j}' - Z'*B{i,j}');
                    cos(theta)*(P'*A{i,j}' + Z'*B{i,j}' - A{i,j}*P - B{i,j}*Z),sin(theta)*(A{i,j}*P + B{i,j}*Z + P*A{i,j}' + Z'*B{i,j}')] <= 0;
                    A{i,j}*P + P*A{i,j}' + B{i,j}*Z + Z'*B{i,j}' - 2*h1*P <=0;
                    A{i,j}*P + P*A{i,j}' + B{i,j}*Z + Z'*B{i,j}' - 2*h2*P >= 0];

        LMI_Hinf{i,j} = [P*A{i,j}' + Z'*B{i,j}' + A{i,j}*P + B{i,j}*Z, Bw{i,j}, P*C{i,j}' + Z'*D{i,j}';
             Bw{i,j}', -delta*eye(size(Bw{i,j},2),size(Bw{i,j},2)), D{i,j}';
             C{i,j}*P + D{i,j}*Z, D{i,j}, -eye(size(D{i,j},1),size(D{i,j},1))] <= 0;

        LMI_H2{i,j} = [[P*A{i,j}' + Z'*B{i,j}' + A{i,j}*P + B{i,j}*Z, P*C{i,j}' + Z'*D{i,j}';
           C{i,j}*P + D{i,j}*Z, -eye(size(C{i,j},1))] <= 0;
                     [J Bw{i,j}';
                        Bw{i,j} P] >= 0;];

        LMI_H2Hinf{i,j} = [LMI_Hinf{i,j};LMI_H2{i,j}];

    end
end

LMI_H2Hinf = [LMI_Hinf;LMI_H2];

% Resultant LMI
LMI = [LMIJ;LMI0;LMI1{1,1};LMI1{1,2};LMI1{2,1};LMI1{2,2};LMI_H2Hinf{1,1};LMI_H2Hinf{1,2};LMI_H2Hinf{2,1};LMI_H2Hinf{2,2}];
% LMI = [LMIJ;LMI0;LMI_H2Hinf{1,1};LMI_H2Hinf{1,2};LMI_H2Hinf{2,1};LMI_H2Hinf{2,2}];

optimize(LMI,trace(J));
checkset(LMI);
Zo = value(Z);
Po = value(P);

% Controller Gain
K = Zo/Po;

%% Results

Gcl = cell(n); % Closed Loop Dynamic Matrices
M = Gcl; % output
Mu = Gcl; % control signal
poles = Gcl; % closed loop poles

for i = 1:n
    for j = 1:n

        % Closed Loop Matrix and poles
        Gcl{i,j} = A{i,j}+B{i,j}*K;
        poles{i,j} = eig(Gcl{i,j});
        M{i,j} = ss((A{i,j}+B{i,j}*K),B{i,j},C{i,j}+D{i,j}*K,D{i,j});
        Mu{i,j} = ss((A{i,j}+B{i,j}*K),B{i,j},K,0);

        figure

        subplot(3,2,1)
        impulse(ss(A{i,j},B{i,j},C{i,j},D{i,j}))
        legend('Open Loop Output')
        title("")
        xlabel("Time","FontSize",8)
        ylabel('y_{OL}(t)')
        set(findall(gcf,'type','line'),'linewidth',3);
        grid on

        subplot(3,2,3)
        impulse(M{i,j});
        legend('Closed Loop Output')
        title("")
        xlabel("Time","FontSize",8)
        ylabel('y_{CL}(t)')
        set(findall(gcf,'type','line'),'linewidth',3);
        grid on

        subplot(3,2,5)
        impulse(Mu{i,j})
        legend('Control Signal')
        title("")
        xlabel("Time","FontSize",8)
        ylabel('u(t)')
        set(findall(gcf,'type','line'),'linewidth',3);
        grid on
        sgtitle("Polytope vertex (" + i +","+ j + ") with $M_P$ =  " + Mp_i{i,j} + " and $K_P$ = " + Kp_i{i,j});
        
        % Plots the pzmap with constraints
        subplot(3,2,[2 4 6])
        t = -tmax:0.1:tmax;
        l1 = -tan(theta)*(-tmax:0.1:0);
        l2 = tan(theta)*(-tmax:0.1:0);
        h_1=h1*ones(size(t));
        h_2=h2*ones(size(t));
        plot(-tmax:0.1:0,l1,'r','LineWidth', 2);hold on; 
        plot(-tmax:0.1:0,l2,'r','LineWidth',2);hold on;
        plot(h_1,t,'r','LineWidth',2); hold on;
        plot(h_2,t,'r','LineWidth',2); hold on;
        scatter(real(poles{i,j}),imag(poles{i,j}),'rx','LineWidth',2);
        hold on
        scatter(real(poles_OL{i,j}),imag(poles_OL{i,j}),'bx','LineWidth',2);
        grid on; 
    end
end

%% Result for a random point in the polytope
% Mp_test = Mp(randi([1 length(Mp)]));
% Kp_test = Kp(randi([1 length(Kp)]));
Mp_test = 5;
Kp_test = 10;


A_test = [0 1 0 0;
    -Kf/Mf -Cf/Mf Kf/Mf Cf/Mf;
    0 0 0 1;
    Kf/Mp_test Cf/Mp_test -(Kp_test+Kf)/Mp_test -(Cp+Cf)/Mp_test];
B_test = B{1,1};
Cz_test = C{1,1};
D_test = D{1,1};

poles_OL_test = eig(A_test);

% Closed Loop Matrix and poles
Gcl_test = A_test+B_test*K;
poles_test = eig(Gcl_test);
M_test = ss((A_test+B_test*K),B_test,Cz_test+D_test*K,D_test);
Mu_test = ss((A_test+B_test*K),B_test,K,0);

figure

subplot(3,2,1)
impulse(ss(A_test,B_test,Cz_test,D_test))
legend('Open Loop Output')
title("")
xlabel("Time","FontSize",8)
ylabel('y_{OL}(t)')
set(findall(gcf,'type','line'),'linewidth',3);
grid on

subplot(3,2,3)
impulse(M_test);
legend('Closed Loop Output')
title("")
xlabel("Time","FontSize",8)
ylabel('y_{CL}(t)')
set(findall(gcf,'type','line'),'linewidth',3);
grid on

subplot(3,2,5)
impulse(Mu_test)
legend('Control Signal')
title("")
xlabel("Time","FontSize",8)
ylabel('u(t)')
set(findall(gcf,'type','line'),'linewidth',3);
grid on

% Plots the pzmap with constraints
subplot(3,2,[2 4 6])
t = -tmax:0.1:tmax;
l1 = -tan(theta)*(-tmax:0.1:0);
l2 = tan(theta)*(-tmax:0.1:0);
h_1=h1*ones(size(t));
h_2=h2*ones(size(t));
plot(-tmax:0.1:0,l1,'r','LineWidth',2);hold on; 
plot(-tmax:0.1:0,l2,'r','LineWidth',2);hold on;
plot(h_1,t,'r','LineWidth',2); hold on;
plot(h_2,t,'r','LineWidth',2); hold on;
scatter(real(poles_test),imag(poles_test),'rx','LineWidth',2);
hold on
scatter(real(poles_OL_test),imag(poles_OL_test),'bx','LineWidth',2);
grid on;

sgtitle("Polytope vertex with $M_p$ =  " + Mp_test + " and $K_p$ = " + Kp_test);

%% H2Hinf plots

w = 1:0.01:10^3;

% Tzw - Open Loop
Tzw = ss(A_test,B_test,Cz_test,D_test);
[mag,~,~] = bode(Tzw,w);
Mag1=squeeze(20*log10(mag(1,1,:)));
% Mag2=squeeze(20*log10(mag(2,1,:)));


% Tzw - Closed Loop
Tzw_cl = ss(A_test+B_test*K,B_test,Cz_test+D_test*K,D_test);
[mag,~,wout] = bode(Tzw_cl,w);
Mag1_cl=squeeze(20*log10(mag(1,1,:)));
% Mag2_cl=squeeze(20*log10(mag(2,1,:)));

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

% subplot(2,1,1)
semilogx(wout,Mag1,wout,Mag1_cl,'LineWidth',2);grid on;
axis([10^0 10^3 -150 0])
ylabel("Magnitude(dB)")
str0 = "Polytope vertex with $M_p$ =  " + Mp_test + " and $K_p$ = " + Kp_test;
str1 = "Open Loop: $||T_{zw}(j \cdot \omega)||_{\infty}$ = ";
str2 = "Closed Loop: $||T_{zw}(j \cdot \omega)||_{\infty}$ = ";
str3 = "Open Loop: $||T_{zw}(j \cdot \omega)||_{2}$ = ";
str4 = "Closed Loop: $||T_{zw}(j \cdot \omega)||_{2}$ = ";
title({str0 ...
    str1 + "$" + Hinf + "$" ...
    str2 + Hinf_cl ...
    str3 + "$" + H2 + "$" ...
    str4 + H2_cl},'Interpreter','latex')

% subplot(2,1,2)
% semilogx(wout,Mag2,wout,Mag2_cl,'LineWidth',2);grid on;
% axis([10^0 10^3 -150 0])
% ylabel("Magnitude(dB)")

legend("$||T_{zw}(j \cdot \omega)||$ - Open Loop","$||T_{zw}(j \cdot \omega)||$ - Closed Loop",'FontSize',10,'Interpreter','latex','Location','best')
xlabel("Frequency $\omega$ (rad/s)",'Interpreter','latex')
