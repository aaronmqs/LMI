%% code
% Author: Aaron Marques
% Date: 05/04/2022
% LMI Control class (Prof. Fabrício Nogueira)

%% Initial configurations
clear; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

%% Process Model
n = 2; % number of uncertainties
c = (450*rand(1,50) + 850); % random numbers
k = (9.96*rand(1,50) + 10); % random numbers
m1 = 288.9;
m2 = 28.58;
kt = 155.900;

H = cell(n); % Polytope Vertices
c_i = H; % uncertainties in the damping coefficient
k_i = H; % uncertainties in the static gain
poles_OL = H; % Open Loop Poles
% State Space Matrices
A = H;
B = H;
C = H;
D = H;

% Polytop plotting
figure
scatter(c,k,'filled');
ylabel("Spring Constant - K")
xlabel("c")
title("Polytopic uncertainties")
hold on
grid on

for i = 1:n
    for j = 1:n
        if mod(j,2) == 0
            c_i{i,j} = max(c);
        else 
            c_i{i,j} = min(c);
        end
        if mod(i,2) == 0
            k_i{i,j} = min(k);
        else 
            k_i{i,j} = max(k);
        end
        scatter(c_i{i,j},k_i{i,j},'rx','LineWidth',2);
        hold on
        A{i,j} = [  0 0 -1 1 0;
                    0 0 0 -1 0;
                    k_i{i,j}/m1 0 -c_i{i,j}/m1 c_i{i,j}/m1 0;
                    -k_i{i,j}/m2 kt/m2 c_i{i,j}/m2 -c_i{i,j}/m2 0;
                    1 0 0 0 0];
        B{i,j} = [0 0 1/m1 -1/m2 0]';
        C{i,j} = [1 0 0 0 0];
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
        step(ss(A{i,j},B{i,j},C{i,j},D{i,j}))
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
        sgtitle("Polytope vertex (" + i +","+ j + ") with c =  " + c_i{i,j} + " and k = " + k_i{i,j});
        
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
k_test = k(randi([1 length(k)]));
c_test = c(randi([1 length(c)]));

A_test = [  0 0 -1 1 0;
            0 0 0 -1 0;
            k_test/m1 0 -c_test/m1 c_test/m1 0;
            -k_test/m2 kt/m2 c_test/m2 -c_test/m2 0;
            1 0 0 0 0];
B_test = [0 0 1/m1 -1/m2 0]';
C_test = [1 0 0 0 0];
D_test = 0;

poles_OL_test = eig(A_test);

% Closed Loop Matrix and poles
Gcl_test = A_test+B_test*K;
poles_test = eig(Gcl_test);
M_test = ss((A_test+B_test*K),B_test,C_test+D_test*K,D_test);
Mu_test = ss((A_test+B_test*K),B_test,K,0);

figure

subplot(3,2,1)
step(ss(A_test,B_test,C_test,D_test))
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
sgtitle("Polytope vertex with $c$ =  " + c_test + " and k = " + k_test);

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
