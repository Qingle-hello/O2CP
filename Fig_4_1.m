clc; clear; close all;

a1 = 0.97822;
a2 = -0.46300;
b1 = log(0.56380);
c2 = -0.00047;

eb1 = exp(b1);

B0 = @(s) (a1 + a2 .* s) ./ (1 + eb1 .* s);
B_m1 = @(s) ((1 - a1) + c2 .* s) ./ (1 + eb1 .* s);

r = @(s) (24 - 6 .* s) ./ (24 + 18 .* s + 6 .* s.^2 + s.^3);

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

denom1 = (1 - abs_alpha) .* (1 - abs_beta);

Nc = 1000;
denom2 = zeros(size(s));
for i = 1:length(s)
    a = alpha(i);
    b = beta(i);
    diff_ab = abs(a - b);
    if diff_ab == 0
        denom2(i) = 0;
    else
        sum_k = 0;
        for k = 0:Nc
            sum_k = sum_k + abs(a^(k+1) - b^(k+1)) / diff_ab;
        end
        denom2(i) = sum_k;
    end
end

figure('Units', 'pixels', 'Position', [100, 100, 600, 600]);
plot(s, abs_alpha, 'LineWidth', 2);
hold on;
plot(s, abs_beta, '-.', 'LineWidth', 2, 'Color', [0.8, 0.4, 0]);
hold off;
set(gca, 'XScale', 'log');
xlim([1e-4, 1e4]);
ylim([0, 1]);
xlabel('$s$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$|\rho_{i}(s)|$', 'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'XTick', [1e-4, 1e-2, 1e0, 1e2, 1e4]);
set(gca, 'FontSize', 30);
set(gca, 'FontName', 'Times New Roman');
legend('$|\rho_1(s)|$', '$|\rho_2(s)|$', 'Interpreter', 'latex', 'FontSize', 30, 'Location', 'best');
grid off;

figure('Units', 'pixels', 'Position', [700, 100, 600, 600]);
hold on;
J_values = [10, 20, 40, 80];
for i = 1:length(J_values)
    J = J_values(i);
    exp_neg_s = r(2 * s ./ J).^(J / 2);
    exp_neg_2s = r(2 * s ./ J).^(J);
    err_num = abs(exp_neg_2s - B0s .* exp_neg_s - B_m1s);
    obj = err_num ./ denom1;
    plot(s, obj, 'LineWidth', 2);
end

exp_neg_s = exp(-s);
exp_neg_2s = exp(-2 .* s);
err_num = abs(exp_neg_2s - B0s .* exp_neg_s - B_m1s);
obj_exp = err_num ./ denom1;
plot(s, obj_exp, '-.', 'LineWidth', 4, 'Color', [0.8235 0.7059 0.5490]);

yline(0.0064, '--', '0.0064', 'FontSize', 20, 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'LineWidth', 2, ...
'FontName', 'Times New Roman');
hold off;
set(gca, 'XScale', 'log');
xlim([1e-4, 1e4]);
ylim([0, 0.02]);
set(gca, 'XTick', [1e-4, 1e-2, 1e0, 1e2, 1e4]);
xlabel('$s$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\gamma_{c}(r,R_1,R_2,J,s)$', 'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 30);
set(gca, 'FontName', 'Times New Roman');
legend_labels = arrayfun(@(x) sprintf('$J=%d$', x), J_values, 'UniformOutput', false);
legend_labels{end+1} = '$\gamma_{e}$';
legend(legend_labels, 'Interpreter', 'latex', 'FontSize', 30, 'Location', 'northwest');
grid off;
box on;

figure('Units', 'pixels', 'Position', [1300, 100, 600, 600]);
hold on;
J_values = [10, 20, 40, 80];
for i = 1:length(J_values)
    J = J_values(i);
    exp_neg_s = r(2 * s ./ J).^(J / 2);
    exp_neg_2s = r(2 * s ./ J).^(J);
    err_num = abs(exp_neg_2s - B0s .* exp_neg_s - B_m1s);
    obj = err_num .* denom2;
    plot(s, obj, 'LineWidth', 2);
end

exp_neg_s = exp(-s);
exp_neg_2s = exp(-2 .* s);
err_num = abs(exp_neg_2s - B0s .* exp_neg_s - B_m1s);
obj_exp = err_num .* denom2;
plot(s, obj_exp, '-.', 'LineWidth', 4, 'Color', [0.8235 0.7059 0.5490]);

yline(0.0062, '--', '0.0062', 'FontSize', 20, 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'LineWidth', 2, ...
'FontName', 'Times New Roman');
hold off;
set(gca, 'XScale', 'log');
xlim([1e-4, 1e4]);
ylim([0, 0.02]);
set(gca, 'XTick', [1e-4, 1e-2, 1e0, 1e2, 1e4]);
xlabel('$s$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\kappa_{c}(r,R_1,R_2,J,s,10^3)$', 'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'FontSize', 30);
set(gca, 'FontName', 'Times New Roman');
legend_labels = arrayfun(@(x) sprintf('$J=%d$', x), J_values, 'UniformOutput', false);
legend_labels{end+1} = '$\kappa_{e}$';
legend(legend_labels, 'Interpreter', 'latex', 'FontSize', 30, 'Location', 'northwest');
grid off;
box on;
