clc; clear; close all;

% BDF2 
a1 = 0.97822;
a2 = -0.46300;
b1 = log(0.56380);
c2 = -0.00047;

eb1 = exp(b1);

B0 = @(s) (a1 + a2 .* s) ./ (1 + eb1 .* s);
B_m1 = @(s) ((1 - a1) + c2 .* s) ./ (1 + eb1 .* s);
% RIIA2
r1 = @(s) (1 - 1/3 .* s)./ (1 + 2/3 .* s + 1/6 .* s.^2);
% RIIA3
r2 = @(s) ((1 - 2/5 .* s + 1/20 .* s.^2) ./ (1 + 3/5 .* s + 3/20 .* s.^2 + 1/60 .* s.^3));


s = logspace(-4, 4, 1000);

B0s = B0(s);
B_m1s = B_m1(s);

sum_ab = B0s;
prod_ab = B_m1s;

disc = sum_ab.^2 + 4 .* prod_ab;

sqrt_disc = sqrt(disc); 

alpha = (sum_ab + sqrt_disc) / 2;
beta = (sum_ab - sqrt_disc) / 2;

abs_alpha = abs(alpha);
abs_beta = abs(beta);


denom = (1 - abs_alpha) .* (1 - abs_beta);


figure('Units', 'pixels', 'Position', [100, 100, 600, 600]);
hold on;
J_values = [10, 20, 40, 80];
max_obj = zeros(length(J_values), 1);
for i = 1:length(J_values)
    J = J_values(i);
    exp_neg_s = r1(2 * s ./ J).^(J / 2);
    exp_neg_2s = r1(2 * s ./ J).^(J);
    err_num = abs(exp_neg_2s - B0s .* exp_neg_s - B_m1s);
    obj = err_num ./ denom;
    plot(s, obj, 'LineWidth', 2);
    max_obj(i) = max(obj);
end


exp_neg_s = exp(-s);
exp_neg_2s = exp(-2 .* s);
err_num = abs(exp_neg_2s - B0s .* exp_neg_s - B_m1s);
obj_exp = err_num ./ denom;
plot(s, obj_exp, 'LineStyle', '-.', 'LineWidth', 4,'Color', [0.8235 0.7059 0.5490]);
max_obj_exp = max(obj_exp);


yline(0.0064, '--', '0.0064', 'FontSize', 20 , 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top','LineWidth', 2,...
'FontName','Times New Roman');
hold off;
set(gca, 'XScale', 'log','FontName','Times New Roman');
xlim([1e-4, 1e4]);
ylim([0, 0.02]);
set(gca, 'XTick', [1e-4, 1e-2, 1e0, 1e2, 1e4]);
xlabel('$s$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\gamma_{c}(r,R_1,R_2,J,s)$', 'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 30);
legend_labels = arrayfun(@(x) sprintf('$J=%d$', x), J_values, 'UniformOutput', false);
legend_labels{end+1} = '$\gamma_{e} (R_1,R_2,s)$';
legend(legend_labels, 'Interpreter', 'latex', 'FontSize', 25, 'Location', 'best');
grid off;
box on;


figure('Units', 'pixels', 'Position', [700, 100, 600, 600]);
hold on;
J_values = [10, 20, 40, 80];
max_obj = zeros(length(J_values), 1);
for i = 1:length(J_values)
    J = J_values(i);
    exp_neg_s = r2(2 * s ./ J).^(J / 2);
    exp_neg_2s = r2(2 * s ./ J).^(J);
    err_num = abs(exp_neg_2s - B0s .* exp_neg_s - B_m1s);
    obj = err_num ./ denom;
    plot(s, obj, 'LineWidth', 2);
    max_obj(i) = max(obj);
end

exp_neg_s = exp(-s);
exp_neg_2s = exp(-2 .* s);
err_num = abs(exp_neg_2s - B0s .* exp_neg_s - B_m1s);
obj_exp = err_num ./ denom;
plot(s, obj_exp, 'LineStyle', '-.', 'LineWidth', 4,'Color', [0.8235 0.7059 0.5490]);
max_obj_exp = max(obj_exp);


yline(0.0064, '--', '0.0064', 'FontSize', 20, 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top','LineWidth', 2,...
'FontName','Times New Roman');
hold off;
set(gca, 'XScale', 'log');
xlim([1e-4, 1e4]);
ylim([0, 0.02]);
set(gca, 'XTick', [1e-4, 1e-2, 1e0, 1e2, 1e4],'FontName','Times New Roman');
xlabel('$s$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\gamma_{c}(r,R_1,R_2,J,s)$', 'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 30);
legend_labels = arrayfun(@(x) sprintf('$J=%d$', x), J_values, 'UniformOutput', false);
legend_labels{end+1} = '$\gamma_{e} (R_1,R_2,s)$';
legend(legend_labels, 'Interpreter', 'latex', 'FontSize',25, 'Location', 'best');
grid off;
box on;


figure('Units', 'pixels', 'Position', [1300, 100, 600, 600]);
J = [2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80];
load('RIIA2.mat');
error_record_RIIA2 = abs(error_record); % Ensure positive for plotting
load('LIIIC3.mat');
error_record_LIIIC3 = abs(error_record);
load('RIIA3.mat');
error_record_RIIA3 = abs(error_record);
loglog(J, error_record_RIIA2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
loglog(J, error_record_LIIIC3, '-^', 'LineWidth', 2, 'MarkerSize', 8);
J1 = [2, 4, 6, 8, 10, 20, 30, 40, 50, 60];
loglog(J1, error_record_RIIA3(1:end-2), '-s', 'LineWidth', 2, 'MarkerSize', 8);

k3 = 3;
c3 = error_record_RIIA2(end) * J(end)^k3;
ref3 = 3*c3 ./ J.^k3;
loglog(J, ref3, '--', 'LineWidth', 2);

k4 = 4;
c4 = error_record_LIIIC3(end) * J(end)^k4;
ref4 = 3*c4 ./ J.^k4;
loglog(J, ref4, '--', 'LineWidth', 2);

k5 = 5;
c5 = error_record_RIIA3(end) * J(end)^k5;
ref5 = 3*c5 ./ J.^k5;
loglog(J, ref5, '--', 'LineWidth', 2);

xlabel('$J$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$|\gamma^{\ast} (r,R_{1},R_{2},J)-\gamma_{e}^{\ast} (R_{1},R_{2})|$', 'Interpreter', 'latex', 'FontSize', 20);

grid off;
xlim([1.5,100]);
set(gca, 'XTick', [2,4,10,20,40,80], 'FontSize', 30,'FontName','Times New Roman');

legend('RIIA2', 'LIIIC3', 'RIIA3', ...
'$\mathcal{O}(J^{-3})$', '$\mathcal{O}(J^{-4})$', '$\mathcal{O}(J^{-5})$', ...
'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 25);
hold off;
