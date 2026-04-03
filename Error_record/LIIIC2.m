clc; clear; close all;

a1 = 0.97822;
a2 = -0.46300;
b1 = log(0.56380);
c2 = -0.00047;

eb1 = exp(b1);

B0 = @(s) (a1 + a2 .* s) ./ (1 + eb1 .* s);
B_m1 = @(s) ((1 - a1) + c2 .* s) ./ (1 + eb1 .* s);

% r = @(s) (-6.*s + 24) ./ (s.^3 + 6.*s.^2 + 18.*s + 24);
r = @(s) (2) ./ (s.^2 + 2*s + 2);

s = logspace(-2, 2, 1000000);

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

figure('Units', 'pixels', 'Position', [100, 100, 1800, 600]);

subplot(1, 3, 1);
plot(s, abs_alpha, 'LineWidth', 4);
set(gca, 'XScale', 'log');
xlim([1e-4, 1e4]);
ylim([0, 1]);
xlabel('$s$', 'Interpreter', 'latex', 'FontSize', 2);
ylabel('$|\rho_1(s)|$', 'Interpreter', 'latex', 'FontSize', 2);
set(gca, 'XTick', [1e-4, 1e-2, 1e0, 1e2, 1e4]);
set(gca, 'FontSize', 30);
grid off;


subplot(1, 3, 2);
plot(s, abs_beta, 'LineWidth', 4);
set(gca, 'XScale', 'log');
xlim([1e-4, 1e4]);
ylim([0, 1]);
set(gca, 'XTick', [1e-4, 1e-2, 1e0, 1e2, 1e4]);
xlabel('$s$', 'Interpreter', 'latex', 'FontSize', 2);
ylabel('$|\rho_2(s)|$', 'Interpreter', 'latex', 'FontSize', 2);
set(gca, 'FontSize', 30);
grid off;


subplot(1, 3, 3);
hold on;
J_values = [2,4,6,8,10,20,30,40,50,60,70,80];

max_obj = zeros(length(J_values), 1);
for i = 1:length(J_values)
    J = J_values(i);
    exp_neg_s = r(2 * s ./ J).^(J / 2);
    exp_neg_2s = r(2 * s ./ J).^(J);
    err_num = abs(exp_neg_2s - B0s .* exp_neg_s - B_m1s);
    denom = (1 - abs_alpha) .* (1 - abs_beta);
    obj = err_num ./ denom;
    plot(s, obj, 'LineWidth', 4);
    max_obj(i) = max(obj);
end


exp_neg_s = exp(-s);
exp_neg_2s = exp(-2 .* s);
err_num = abs(exp_neg_2s - B0s .* exp_neg_s - B_m1s);
denom = (1 - abs_alpha) .* (1 - abs_beta);
obj_exp = err_num ./ denom;
plot(s, obj_exp, 'LineStyle', '-.', 'LineWidth', 6,'Color', [0.8235 0.7059 0.5490]);
max_obj_exp = max(obj_exp);


yline(0.0064, '--', '0.0064', 'FontSize', 20, 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top','LineWidth', 2);
hold off;
set(gca, 'XScale', 'log');
xlim([1e-4, 1e4]);
ylim([0, 0.03]);
set(gca, 'XTick', [1e-4, 1e-2, 1e0, 1e2, 1e4]);
xlabel('$s$', 'Interpreter', 'latex', 'FontSize', 2);
ylabel('$\gamma_{c}(r,R_1,R_2,J,s)$', 'Interpreter', 'latex', 'FontSize', 2);
set(gca, 'FontSize', 30);
legend_labels = arrayfun(@(x) sprintf('$J=%d$', x), J_values, 'UniformOutput', false);
legend_labels{end+1} = '$\gamma_{e} (R_1,R_2,s)$';
legend(legend_labels, 'Interpreter', 'latex', 'FontSize', 30, 'Location', 'best');
grid off;
box on;

error_record = [abs(max_obj - max_obj_exp)];
save('LIIIC2.mat', 'error_record');
