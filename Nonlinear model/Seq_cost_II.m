% We test for the sequential running time
% SDIRK2, BDF2, OCP, TS OCP.
clc;clear;close all;
addpath('Solver');

% J=20/2;cL=1;
% J=50/2;cL=1;
% J=20/2;cL=5;
% J=50/2;cL=5;

% J=20/2;cL=10;
J=50/2;cL=10;

CP1 = @SDIRK2;
CP2 = @BDF2_solver;
CP3 = @OCP;
CP4 = @MOCP1_solver;
CP6 = @MOCP1_solver_Extplt_nl;

FP = @LIIIC3_solver;

% Parareal setting
Tend = 10; tau = 0.01/cL;
Tau = tau*J;Nc = Tend/Tau;
iter_max = 20;

%% Preparation
N = 1000-1; dx = 1 / (N + 1);
x = (1:N)' * dx; % Grid points

A = spdiags([-ones(N,1), 2*ones(N,1), -ones(N,1)], -1:1, N, N) / dx^2;

u_exact = @(x,t) sin(pi * x) .* cos(pi * t);
g = @(x,t) -pi * sin(pi * x) .* sin(pi * t) ...
          + pi^2 * sin(pi * x) .* cos(pi * t) ...
          - (cL * u_exact(x,t) .* (1 - u_exact(x,t).^2));

% f(u,t) = cL * u.^3 + g(x,t), dfdu(u,t) = cL * diag(3 u.^2)

f_func = @(u,t) cL * u.*(1-u.^2) + g(x,t);
dfdu_func = @(u,t) spdiags(cL * (1- 3*u.^2), 0, N, N);

u0 = u_exact(x, 0);

u_seq=zeros(N,Nc+1);u_seq(:,1)=u0;
u_para=zeros(N,Nc+1);u_para(:,1)=u0;
u_para_new=zeros(N,Nc+1);u_para_new(:,1)=u0;
u_temp=zeros(N,Nc+1);
err=zeros(iter_max+1,1);

%% compute the first value
for i = 1:2
    tspan = [(i-1)*Tau i*Tau];
    [~, u_new] = FP(A, f_func, dfdu_func, u_seq(:,i), tspan, tau);
    u_seq(:,i+1) = u_new(:,end);
end


%% BDF2
times1 = zeros(50,1);
for run=1:50
    u_para1 = u_para;
    u_para1(:,2)=u_seq(:,2);
    t1 = tic;
    % Initialization
    for i=2:Nc
        tspan = [(i-1)*Tau i*Tau];
        [~,u_new] = CP2(A, f_func, dfdu_func, [u_para(:,i-1),u_para(:,i)], tspan, Tau);
        u_para1(:,i+1) = u_new(:,end);
    end
    times1(run) = toc(t1);
end
avg_time1 = mean(times1);

fprintf('BDF2:\n Average Time is %.4f.\n ',avg_time1);

%% O2CP
times3 = zeros(50,1);
for run=1:50
    u_para3 = u_para;
    u_para3(:,2)=u_seq(:,2);
    t1 = tic;

    for i=2:Nc
        tspan = [(i-1)*Tau i*Tau];
        [~,u_new] = CP4(A, f_func, dfdu_func, [u_para(:,i-1),u_para(:,i)], tspan, Tau);
        u_para3(:,i+1) = u_new(:,end);
    end
    times3(run) = toc(t1);
end
avg_time3 = mean(times3);

fprintf('O2CP:\n Average Time is %.4f.\n',avg_time3);


%% O2CP-E,
% The function: MOCP1_solver_Extplt.m
% and the function: MOCP1_solver_Extplt_nl.m,
% the first uses 'fsolve' and the second uses 'A\b'.
% We test the time for MOCP1_solver_Extplt_nl.m

times5 = zeros(50,1);
for run=1:50
    u_para5 = u_para;
    u_para5(:,2)=u_seq(:,2);
    t1 = tic;
    % Initialization
    for i=2:Nc
        tspan = [(i-1)*Tau i*Tau];
        [~,u_new] = CP6(A, f_func, dfdu_func, [u_para(:,i-1),u_para(:,i)], tspan, Tau);
        u_para5(:,i+1) = u_new(:,end);
    end
    times5(run) = toc(t1);
end
avg_time5 = mean(times5);

fprintf('O2CP-E:\n Average Time is %.4f.\n',avg_time5);
