%% Example of control for a second-order process with polytopic uncertainties with spacifications defined by D-Stability LMI conditions
% Author: Aaron Marques
% Date: 04/29/2022
% LMI Control class (Prof. Fabrício Nogueira)

%% Initial configurations
clear; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

%% Process Model
n = 2; % number of uncertainties
delta = -2*rand(1,50) + 1; % random numbers in (-1,1): uncertainties in the damping coefficient
k = -0.2*rand(1,50) + 1.1; % random numbers in (0.9,1.1): uncertainties in the static gain

H = cell(n); % Polytope Vertices
delta_i = H; % uncertainties in the damping coefficient
k_i = H; % uncertainties in the static gain
poles_OL = H; % Open Loop Poles
% State Space Matrices
A = H;
B = H;
C = H;
D = H;

% Polytop plotting
figure
scatter(delta,k,'filled');
ylabel("Static gain - K")
xlabel("Damping coefficient depending factor")
title("Polytopic uncertainties")
hold on
grid on

for i = 1:n
    for j = 1:n
        if mod(j,2) == 0
            delta_i{i,j} = max(delta);
        else 
            delta_i{i,j} = min(delta);
        end
        if mod(i,2) == 0
            k_i{i,j} = min(k);
        else 
            k_i{i,j} = max(k);
        end
        scatter(delta_i{i,j},k_i{i,j},'rx','LineWidth',2);
        hold on
        H{i,j} = tf(k_i{i,j},[1 0.1*delta_i{i,j} 2]);
        [nH,dH] = tfdata(H{i,j},'v');
        [nH,dH] = eqtflength(nH,dH);
        [Ai,Bi,Ci,Di] = tf2ss(nH,dH);
        A{i,j} = Ai;
        B{i,j} = Bi;
        C{i,j} = Ci;
        D{i,j} = Di;
        poles_OL{i,j} = eig(Ai);
    end
end

%% Performance parameters
t_min = 2; % min settling time
t_max = 5; % max settling time
UP = 0.1; % overshoot
zeta = sqrt((log(UP))^2/(pi^2 + (log(UP))^2)); % damping coefficient
h1 = -4/t_max; % Vertical Strip
h2 = -4/t_min; % Vertical Strip
theta = acos(zeta); % Conic Sector angle

%% LMI (Controller Synthesis)

% Decision Variables
P = sdpvar(size(A{1,1},1),size(A{1,1},1));
Z = sdpvar(size(B{1,1},2),size(A{1,1},1));

LMI0 = P>=0;

LMI1 = cell(n); % Matrix of LMIs

% Conic Sector Constraints
% Left Vertical Strip Constraint
% Right Vertical Strip Constraint

for i = 1:n
    for j = 1:n
        LMI1{i,j} = [[sin(theta)*(A{i,j}*P + B{i,j}*Z + P*A{i,j}' + Z'*B{i,j}'),cos(theta)*(A{i,j}*P + B{i,j}*Z - P*A{i,j}' - Z'*B{i,j}');
                    cos(theta)*(P'*A{i,j}' + Z'*B{i,j}' - A{i,j}*P - B{i,j}*Z),sin(theta)*(A{i,j}*P + B{i,j}*Z + P*A{i,j}' + Z'*B{i,j}')] <= 0;
                    A{i,j}*P + P*A{i,j}' + B{i,j}*Z + Z'*B{i,j}' - 2*h1*P <=0;
                    A{i,j}*P + P*A{i,j}' + B{i,j}*Z + Z'*B{i,j}' - 2*h2*P >= 0];
    end
end

% Resultant LMI
LMI = [LMI0;LMI1{1,1};LMI1{1,2};LMI1{2,1};LMI1{2,2}];

%ops = sdpsettings('solver','gurobi');
optimize(LMI);
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
        step(H{i,j})
        legend('Open Loop Output')
        title("")
        xlabel("Time","FontSize",8)
        ylabel('y_{OL}(t)')
        set(findall(gcf,'type','line'),'linewidth',3);
        grid on

        subplot(3,2,3)
        step(M{i,j});
        legend('Closed Loop Output')
        title("")
        xlabel("Time","FontSize",8)
        ylabel('y_{CL}(t)')
        set(findall(gcf,'type','line'),'linewidth',3);
        grid on

        subplot(3,2,5)
        step(Mu{i,j})
        legend('Control Signal')
        title("")
        xlabel("Time","FontSize",8)
        ylabel('u(t)')
        set(findall(gcf,'type','line'),'linewidth',3);
        grid on
        sgtitle("Polytope vertix with $\delta$ =  " + delta_i{i,j} + " and k = " + k_i{i,j});
        
        % Plots the pzmap with constraints
        subplot(3,2,[2 4 6])
        t = -200:0.1:200;
        l1 = -theta*t;
        l2 = theta*t;
        h_1=h1*ones(size(t));
        h_2=h2*ones(size(t));
        plot(t,l1,'r','LineWidth', 2);
        hold on; plot(t,l2,'r','LineWidth',2);
        hold on;plot(h_1,t,'r','LineWidth',2); hold on;
        hold on;plot(h_2,t,'r','LineWidth',2); hold on;
        scatter(real(poles{i,j}),imag(poles{i,j}),'rx','LineWidth',2);
        hold on
        scatter(real(poles_OL{i,j}),imag(poles_OL{i,j}),'bx','LineWidth',2);
        grid on; 
        % Configures the limits of the axis for better visual
        xmin = 1.5*h2;
        x_aux = real(poles_OL{i,j}) + 0.1*abs(real(poles_OL{i,j}));
        xmax = x_aux(1);
        y_aux = max(abs(imag(poles{i,j})),abs(imag(poles_OL{i,j})));
        ymin = -y_aux(1);
        ymax = y_aux(1);
        axis([xmin xmax ymin ymax])
    end
end

%% Result for a random point in the polytope
k_test = k(randi([1 length(k)]));
delta_test = delta(randi([1 length(delta)]));
H_test = tf(k_test,[1 0.1*delta_test 2]); % test matrix
[nH,dH] = tfdata(H_test,'v');
[nH,dH] = eqtflength(nH,dH);
[A_test,B_test,C_test,D_test] = tf2ss(nH,dH);
poles_OL_test = eig(A_test);

% Closed Loop Matrix and poles
Gcl_test = A_test+B_test*K;
poles_test = eig(Gcl_test);
M_test = ss((A_test+B_test*K),B_test,C_test+D_test*K,D_test);
Mu_test = ss((A_test+B_test*K),B_test,K,0);

figure

subplot(3,2,1)
step(H_test)
legend('Open Loop Output')
title("")
xlabel("Time","FontSize",8)
ylabel('y_{OL}(t)')
set(findall(gcf,'type','line'),'linewidth',3);
grid on

subplot(3,2,3)
step(M_test);
legend('Closed Loop Output')
title("")
xlabel("Time","FontSize",8)
ylabel('y_{CL}(t)')
set(findall(gcf,'type','line'),'linewidth',3);
grid on

subplot(3,2,5)
step(Mu_test)
legend('Control Signal')
title("")
xlabel("Time","FontSize",8)
ylabel('u(t)')
set(findall(gcf,'type','line'),'linewidth',3);
grid on
sgtitle("Polytope vertix with $\delta$ =  " + delta_test + " and k = " + k_test);

% Plots the pzmap with constraints
subplot(3,2,[2 4 6])
t = -200:0.1:200;
l1 = -tan(theta)*t;
l2 = tan(theta)*t;
h_1=h1*ones(size(t));
h_2=h2*ones(size(t));
plot(t,l1,'r','LineWidth', 2);
hold on; plot(t,l2,'r','LineWidth',2);
hold on;plot(h_1,t,'r','LineWidth',2); hold on;
hold on;plot(h_2,t,'r','LineWidth',2); hold on;
scatter(real(poles_test),imag(poles_test),'rx','LineWidth',2);
hold on
scatter(real(poles_OL_test),imag(poles_OL_test),'bx','LineWidth',2);
grid on; 
% Configures the limits of the axis for better visual
xmin = 1.5*h2;
x_aux = real(poles_OL_test) + 0.1*abs(real(poles_OL_test));
xmax = x_aux(1);
y_aux = max(abs(imag(poles_test)),abs(imag(poles_OL_test)));
ymin = -y_aux(1);
ymax = y_aux(1);
axis([xmin xmax ymin ymax])

