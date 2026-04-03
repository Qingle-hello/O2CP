% We test for the sequential running time
% In the parareal test, we use the general solvers in 'Solver' (fsolve).
% for linear problems, we use the linear solvers in 'Solver_linear' (A\b).

clc;clear;close all;
addpath('Solver_linear');
J=50;
CP1 = @SDIRK2_l;
CP2 = @BDF2_solver_l;
CP3 = @OCP;
CP4 = @MOCP1_solver_l;
FP = @LIIIC3_solver;

% Parareal setting
Tend = 10; tau = 0.01;
% J=5;
Tau = tau*J;Nc = Tend/Tau;
iter_max = 20;

%% Preparation
N = 1000-1; dx = 1 / (N + 1);
x = (1:N)' * dx; % Grid points
% A = -Delta^h ≈ tridiag(-1,2,-1)/dx^2
A = spdiags([-ones(N,1), 2*ones(N,1), -ones(N,1)], -1:1, N, N) / dx^2;
% Exact solution and g
u_exact = @(x,t) sin(pi * x) .* cos(pi * t);
g = @(x,t) -pi * sin(pi * x) .* sin(pi * t) ...
          + pi^2 * sin(pi * x) .* cos(pi * t);


f_func = @(u,t) 0.*g(x,t);
dfdu_func = @(u,t) spdiags(0 * u.^2, 0, N, N);

u0 = 1.*(x<0.5);
u_seq=zeros(N,Nc+1);u_seq(:,1)=u0;
u_para=zeros(N,Nc+1);u_para(:,1)=u0;
u_para_new=zeros(N,Nc+1);u_para_new(:,1)=u0;
u_temp=zeros(N,Nc+1);
err=zeros(iter_max+1,1);
%% The reference solution
for i=1:Nc
    tspan = [(i-1)*Tau i*Tau];
    [~,u_new] = FP(A, f_func, dfdu_func, u_seq(:,i), tspan, tau);
    u_seq(:,i+1) = u_new(:,end);
end
fprintf('Reference solution is computed.\n\n Starting Parareal...\n');
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

% Test L^2 error 
diff = u_para1 - u_seq;
err(1) = max(sqrt(dx * sum(diff.^2, 1)));
fprintf('SDIRK2:\n Average Time is %.4f.\n Error is %.4f.\n',avg_time1,err(1));

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
fprintf('OCP:\n Average Time is %.4f.\n Error is %.4f.\n',avg_time3,err(1));
