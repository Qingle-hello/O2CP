clc; clear; close all;
addpath('Solver');

CP = @ LIIIC2_solver;
FP = @ RIIA3_solver;

Tend = 2; tau = 0.01; J=5; Tau = tau*J; Nc = Tend/Tau;
iter_max = 10;

N = 100-1; dx = 1 / (N + 1);
x = (1:N)' * dx;
A = spdiags([-ones(N,1), 2*ones(N,1), -ones(N,1)], -1:1, N, N) / dx^2;

u_exact = @(x,t) sin(pi * x) .* cos(pi * t);
g = @(x,t) -pi * sin(pi * x) .* sin(pi * t) ...
          + pi^2 * sin(pi * x) .* cos(pi * t) ...
          - (sin(pi * x) .* cos(pi * t)).^3;

f_func = @(u,t) u.^3 + g(x,t);
dfdu_func = @(u,t) spdiags(3 * u.^2, 0, N, N);
u0 = sqrt(x+1).*sin(2*pi*x);

u_seq=zeros(N,Nc+1);u_seq(:,1)=u0;
u_para=zeros(N,Nc+1);u_para(:,1)=u0;
u_para_new=zeros(N,Nc+1);u_para_new(:,1)=u0;
u_temp=zeros(N,Nc+1);
err=zeros(iter_max+1,1);

for i=1:Nc
    tspan = [(i-1)*Tau i*Tau];
    [~,u_new] = FP(A, f_func, dfdu_func, u_seq(:,i), tspan, tau);
    u_seq(:,i+1) = u_new(:,end);
end
fprintf('Reference solution is computed.\n\n Starting Parareal...\n');

for i=1:Nc
    tspan = [(i-1)*Tau i*Tau];
    [~,u_new] = CP(A, f_func, dfdu_func, u_para(:,i), tspan, Tau);
    u_para(:,i+1) = u_new(:,end);
end

diff = u_para - u_seq;
err(1) = max(sqrt(dx * sum(diff.^2, 1)));
fprintf('Initialization completed. Error is %.4f.\n',err(1));

total_tic = tic;

for m=1:iter_max
    iter_tic = tic;

    for i=1:Nc
        tspan = [(i-1)*Tau i*Tau];
        [~,u_new] = FP(A, f_func, dfdu_func, u_para(:,i), tspan, tau);
        u_temp(:,i+1) = u_new(:,end);
    end

    for i = 1:Nc
        tspan = [(i-1)*Tau i*Tau];
        [~, uc_new] = CP(A, f_func, dfdu_func, u_para_new(:,i), tspan, Tau);
        [~, uc_old] = CP(A, f_func, dfdu_func, u_para(:,i), tspan, Tau);
        u_para_new(:,i+1) = uc_new(:,end) + u_temp(:,i+1) - uc_old(:,end);
    end

    u_para = u_para_new;

    diff = u_para - u_seq;
    err(m+1) = max(sqrt(dx * sum(diff.^2, 1)));

    iter_time = toc(iter_tic);
    total_time_so_far = toc(total_tic);
    avg_time_per_iter = total_time_so_far / m;
    estimated_remaining = avg_time_per_iter * (iter_max - m);
    fprintf('Iteration %d completed. Estimated time remaining: %.2f seconds.\n',...
        m, estimated_remaining);
end

iter1 = 2; iter2 = 8;
alpha = exp(log(err(iter1+1:iter2+1)) - log(err(iter1:iter2)));
avg_alpha = mean(alpha);
fprintf('\n\nAverage convergence rate from iteration %d to %d: %.4f\n', iter1, iter2-1, avg_alpha);

figure(1);
clf;
m_list = 0:iter_max;
semilogy(m_list, err, '-o', ...
'LineWidth', 1.8, ...
'MarkerSize', 6, ...
'MarkerFaceColor', [0 0.4470 0.7410], ...
'Color', [0 0.4470 0.7410], ...
'DisplayName', 'Parareal Error');
hold on;
rho = 0.08;
ref1 = rho .^ m_list;
max_plot_iter = find(ref1 > min(err)/10, 1, 'last');
semilogy(m_list(1:max_plot_iter), ref1(1:max_plot_iter), '--', ...
'LineWidth', 1.5, ...
'Color', [0.85 0.325 0.098], ...
'DisplayName', sprintf('\\rho = %.2f', rho));
hold off;

xlabel('Iteration Number', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Maximum L^2 Error', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast', ...
'FontSize', 20, ...
'Box', 'on', ...
'Interpreter', 'tex');
title('Convergence History of Parareal Method (LIIIC2)', ...
'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);
set(gca, 'LineWidth', 1.2);
set(gca, 'Box', 'on');
set(gca, 'MinorGridLineStyle', '-');
set(gca, 'MinorGridAlpha', 0.3);

tight_lines = get(gca, 'TightInset');
set(gca, 'Position', [tight_lines(1) tight_lines(2) ...
                      1 - tight_lines(1) - tight_lines(3) - 0.08 ...
                      1 - tight_lines(2) - tight_lines(4)]);

if min(err) < 1e-10
    ytickformat('%.1e');
end

cp_name = func2str(CP);
filename = [cp_name '_J_' num2str(J) '.mat'];
dir_path = 'Fig1_J_4';
save(fullfile(dir_path, filename), 'err');
