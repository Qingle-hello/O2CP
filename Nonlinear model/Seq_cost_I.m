% We test for the sequential running time
% SDIRK2, BDF2, OCP, O2CP, O2CP-E.
clc;clear;close all;
addpath('Solver');
% J=20;cL=1;
% J=50;cL=1;
% J=20;cL=5;
% J=50;cL=5;

% J=20;cL=10;
J=50;cL=10;

CP1 = @SDIRK2;
CP2 = @BDF2_solver;
CP3 = @OCP;
CP4 = @MOCP1_solver;
FP = @LIIIC3_solver;

Tend = 10; tau = 0.01/cL;

Tau = tau*J;Nc = Tend/Tau;
iter_max = 20;

%% Preparation
N = 1000-1; dx = 1 / (N + 1);
x = (1:N)' * dx;

A = spdiags([-ones(N,1), 2*ones(N,1), -ones(N,1)], -1:1, N, N) / dx^2;

u_exact = @(x,t) sin(pi * x) .* cos(pi * t);
g = @(x,t) -pi * sin(pi * x) .* sin(pi * t) ...
          + pi^2 * sin(pi * x) .* cos(pi * t) ...
          - (cL * u_exact(x,t) .* (1 - u_exact(x,t).^2));


f_func = @(u,t) cL * u.*(1-u.^2) + g(x,t);
dfdu_func = @(u,t) spdiags(cL * (1- 3*u.^2), 0, N, N);


u0 = u_exact(x, 0);

u_seq=zeros(N,Nc+1);u_seq(:,1)=u0;
u_para=zeros(N,Nc+1);u_para(:,1)=u0;
u_para_new=zeros(N,Nc+1);u_para_new(:,1)=u0;
u_temp=zeros(N,Nc+1);
err=zeros(iter_max+1,1);



%% SDIRK2
times1 = zeros(50,1);
for run=1:50
    u_para1 = u_para;
    u_para1(:,2)=u_seq(:,2);
    t1 = tic;
    % Initialization
    for i=2:Nc
        tspan = [(i-1)*Tau i*Tau];
        [~,u_new] = CP1(A, f_func, dfdu_func, u_para1(:,i), tspan, Tau);
        u_para1(:,i+1) = u_new(:,end);
    end
    times1(run) = toc(t1);
end
avg_time1 = mean(times1);

fprintf('SDIRK2:\n Average Time is %.4f\n',avg_time1);

%% OCP
times3 = zeros(50,1);
for run=1:50
    u_para3 = u_para;
    u_para3(:,2)=u_seq(:,2);
    t1 = tic;
    % Initialization
    for i=2:Nc
        tspan = [(i-1)*Tau i*Tau];
        [~,u_new] = CP3(A, f_func, dfdu_func, u_para3(:,i), tspan, Tau);
        u_para3(:,i+1) = u_new(:,end);
    end
    times3(run) = toc(t1);
end
avg_time3 = mean(times3);

diff = u_para3 - u_seq;
err(1) = max(sqrt(dx * sum(diff.^2, 1)));
fprintf('OCP:\n Average Time is %.4f.\n',avg_time3);
