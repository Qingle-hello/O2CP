clc; clear; close all;

addpath('Solver_linear');

J = 50;
CP1 = @SDIRK2_l;
CP2 = @BDF2_solver_l;
CP3 = @OCP;
CP4 = @MOCP1_solver_l;

% Parareal setting
Tend = 10; 
tau = 0.01;
Tau = tau * J; 
Nc = Tend / Tau;

%% Preparation
N = 1000 - 1; 
dx = 1 / (N + 1);
x = (1:N)' * dx;

% Laplacian matrix A = -Delta^h
A = spdiags([-ones(N,1), 2*ones(N,1), -ones(N,1)], -1:1, N, N) / dx^2;

% Source term (仅用于构造 f)
g = @(x,t) -pi * sin(pi * x) .* sin(pi * t) + pi^2 * sin(pi * x) .* cos(pi * t);

% 本次测试为线性情况（f 与 u 无关）
f_func   = @(u,t) 0 .* g(x,t);                    % f = 0
dfdu_func = @(u,t) spdiags(zeros(N,1), 0, N, N); % 零 Jacobian（即使不使用也必须传入）

% 初始条件
u0 = 1 .* (x < 0.5);

%% ====================== Reference Solution (LIIIC3 linear) ======================
FP_l = @LIIIC3_solver_l;   % 线性专用版本（不需要 dfdu）

fprintf('Computing reference solution with LIIIC3_solver_l (linear, dfdu=0)...\n');

times_ref = zeros(50, 1);
u_seq = zeros(N, Nc + 1);
u_seq(:,1) = u0;

for run = 1:50
    u_temp = zeros(N, Nc + 1);
    u_temp(:,1) = u0;
    
    t_ref = tic;
    for i = 1:Nc
        tspan = [(i-1)*Tau, i*Tau];
        [~, u_new] = FP_l(A, f_func, u_temp(:,i), tspan, tau);
        u_temp(:,i+1) = u_new(:,end);
    end
    times_ref(run) = toc(t_ref);
    
    if run == 50
        u_seq = u_temp;               % 保留最终参考解
    end
end

avg_time_ref = mean(times_ref);

fprintf('LIIIC3_solver_l (Reference):\n');
fprintf('   Average Time over 50 runs: %.4f seconds\n', avg_time_ref);
fprintf('   Reference solution computed.\n\n Starting Parareal coarse propagator tests...\n');

%% ====================== SDIRK2 ======================
times1 = zeros(50,1);
for run = 1:50
    u_para1 = zeros(N, Nc+1);
    u_para1(:,1) = u0;
    u_para1(:,2) = u_seq(:,2);
    
    t1 = tic;
    for i = 2:Nc
        tspan = [(i-1)*Tau, i*Tau];
        [~, u_new] = CP1(A, f_func, dfdu_func, u_para1(:,i), tspan, Tau);
        u_para1(:,i+1) = u_new(:,end);
    end
    times1(run) = toc(t1);
end
avg_time1 = mean(times1);

diff = u_para1 - u_seq;
err1 = max(sqrt(dx * sum(diff.^2, 1)));

fprintf('SDIRK2:\n');
fprintf('   Average Time: %.4f s\n', avg_time1);
fprintf('   Max L^2 Error: %.4e\n\n', err1);

%% ====================== OCP ======================
times3 = zeros(50,1);
for run = 1:50
    u_para3 = zeros(N, Nc+1);
    u_para3(:,1) = u0;
    u_para3(:,2) = u_seq(:,2);
    
    t3 = tic;
    for i = 2:Nc
        tspan = [(i-1)*Tau, i*Tau];
        [~, u_new] = CP3(A, f_func, dfdu_func, u_para3(:,i), tspan, Tau);
        u_para3(:,i+1) = u_new(:,end);
    end
    times3(run) = toc(t3);
end
avg_time3 = mean(times3);

diff = u_para3 - u_seq;
err3 = max(sqrt(dx * sum(diff.^2, 1)));

fprintf('OCP:\n');
fprintf('   Average Time: %.4f s\n', avg_time3);
fprintf('   Max L^2 Error: %.4e\n\n', err3);

% 如需继续测试 BDF2 / MOCP1，只需复制上面模式（已准备好 CP2、CP4）
