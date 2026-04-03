clc;clear;close all;
J = [2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80];

load('RIIA2.mat');
error_record_RIIA2 = abs(error_record);  % Ensure positive for plotting

load('LIIIC3.mat');
error_record_LIIIC3 = abs(error_record);

load('RIIA3.mat');
error_record_RIIA3 = abs(error_record);

figure('Position', [100, 100, 800, 600]);

loglog(J, error_record_RIIA2, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
loglog(J, error_record_LIIIC3, '-^', 'LineWidth', 1.5, 'MarkerSize', 6);
loglog(J, error_record_RIIA3, '-s', 'LineWidth', 1.5, 'MarkerSize', 6);

k3 = 3;
c3 = error_record_RIIA2(end) * J(end)^k3;
ref3 = 3*c3 ./ J.^k3;
loglog(J, ref3, '--', 'LineWidth', 1);

k4 = 4;
c4 = error_record_LIIIC3(end) * J(end)^k4;
ref4 = 3*c4 ./ J.^k4;
loglog(J, ref4, '--', 'LineWidth', 1);

k5 = 5;
c5 = error_record_RIIA3(end) * J(end)^k5;
ref5 = 3*c5 ./ J.^k5;
loglog(J, ref5, '--', 'LineWidth', 1);

xlabel('$J$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$|\gamma^{\ast} (r,R_{1},R_{2},J)-\gamma_{e}^{\ast} (R_{1},R_{2})|$', 'Interpreter', 'latex', 'FontSize', 14);

legend('RIIA2', 'LIIIC3', 'RIIA3', ...
       '$\mathcal{O}(J^{-3})$', '$\mathcal{O}(J^{-4})$', '$\mathcal{O}(J^{-5})$', ...
       'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 12);


grid off;
xlim([1.5,100]);
set(gca, 'XTick', [2,4,10,20,40,80], 'FontSize', 12);
hold off;
