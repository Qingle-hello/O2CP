clc; clear; close all;

iter_max = 10;
iter = 0:iter_max;

addpath('data_plot/Fig2_nonsmooth');

load('BDF2_solver_J_20.mat', 'err');
err_BDF2_J20_ns = err;
load('MOCP1_solver_J_20.mat', 'err');
err_MOCP1_J20_ns = err;
load('OCP_J_20.mat', 'err');
err_OCP_J20_ns = err;
load('SDIRK2_J_20.mat', 'err');
err_SDIRK2_J20_ns = err;

load('BDF2_solver_J_50.mat', 'err');
err_BDF2_J50_ns = err;
load('MOCP1_solver_J_50.mat', 'err');
err_MOCP1_J50_ns = err;
load('OCP_J_50.mat', 'err');
err_OCP_J50_ns = err;
load('SDIRK2_J_50.mat', 'err');
err_SDIRK2_J50_ns = err;

rmpath('data_plot/Fig2_nonsmooth');

addpath('data_plot/Fig2_smooth');

load('BDF2_solver_J_20.mat', 'err');
err_BDF2_J20_s = err;
load('MOCP1_solver_J_20.mat', 'err');
err_MOCP1_J20_s = err;
load('OCP_J_20.mat', 'err');
err_OCP_J20_s = err;
load('SDIRK2_J_20.mat', 'err');
err_SDIRK2_J20_s = err;

load('BDF2_solver_J_50.mat', 'err');
err_BDF2_J50_s = err;
load('MOCP1_solver_J_50.mat', 'err');
err_MOCP1_J50_s = err;
load('OCP_J_50.mat', 'err');
err_OCP_J50_s = err;
load('SDIRK2_J_50.mat', 'err');
err_SDIRK2_J50_s = err;

rmpath('data_plot/Fig2_smooth');

fig1 = figure('Position', [100, 100, 600, 400]);
hold on;
plot(iter, err_BDF2_J20_ns, '-o', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'BDF2 ');
plot(iter, err_MOCP1_J20_ns, '-s', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'O2CP');
plot(iter, err_OCP_J20_ns, '-^', 'LineWidth',3, 'MarkerSize', 12, 'DisplayName', 'OCP');
plot(iter, err_SDIRK2_J20_ns, '-d', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'SDIRK2');
plot(iter, 0.05 * 0.151.^iter, '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'DisplayName', sprintf('Ref=$0.151$'));
plot(iter, 0.05 * 0.0062.^iter, '--', 'LineWidth', 2, 'Color', [0.6 0.4 0.4], 'DisplayName', sprintf('Ref=$0.0062$'));
ylim([1e-12,2]);
set(gca, 'YScale', 'log');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
xlabel('Iteration $k$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$e$', 'Interpreter', 'latex', 'FontSize', 20);
legend('Interpreter', 'latex', 'FontSize', 16, 'Location', 'northeast');
box on;
grid off;
hold off;

fig2 = figure('Position', [700, 100, 600, 400]);
hold on;
plot(iter, err_BDF2_J50_ns, '-o', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'BDF2 ');
plot(iter, err_MOCP1_J50_ns, '-s', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'O2CP');
plot(iter, err_OCP_J50_ns, '-^', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'OCP');
plot(iter, err_SDIRK2_J50_ns, '-d', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'SDIRK2');
plot(iter, 0.05 * 0.151.^iter, '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'DisplayName', sprintf('Ref=$0.151$'));
plot(iter, 0.05 * 0.0062.^iter, '--', 'LineWidth', 2, 'Color', [0.6 0.4 0.4], 'DisplayName', sprintf('Ref=$0.0062$'));
ylim([1e-12,2]);
set(gca, 'YScale', 'log');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
xlabel('Iteration $k$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$e$', 'Interpreter', 'latex', 'FontSize', 20);
legend('Interpreter', 'latex', 'FontSize', 16, 'Location', 'northeast');
box on;
grid off;
hold off;

fig3 = figure('Position', [100, 500, 600, 400]);
hold on;
plot(iter, err_BDF2_J20_s, '-o', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'BDF2 ');
plot(iter, err_MOCP1_J20_s, '-s', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'O2CP');
plot(iter, err_OCP_J20_s, '-^', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'OCP');
plot(iter, err_SDIRK2_J20_s, '-d', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'SDIRK2');
plot(iter, 0.05 * 0.151.^iter, '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'DisplayName', sprintf('Ref=$0.151$'));
plot(iter, 0.05 * 0.0062.^iter, '--', 'LineWidth', 2, 'Color', [0.6 0.4 0.4], 'DisplayName', sprintf('Ref=$0.0062$'));
ylim([1e-12,2]);
set(gca, 'YScale', 'log');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
xlabel('Iteration $k$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$e$', 'Interpreter', 'latex', 'FontSize', 20);
legend('Interpreter', 'latex', 'FontSize', 16, 'Location', 'northeast');
box on;
grid off;
hold off;

fig4 = figure('Position', [700, 500, 600, 400]);
hold on;
plot(iter, err_BDF2_J50_s, '-o', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'BDF2 ');
plot(iter, err_MOCP1_J50_s, '-s', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'O2CP');
plot(iter, err_OCP_J50_s, '-^', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'OCP');
plot(iter, err_SDIRK2_J50_s, '-d', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'SDIRK2');
plot(iter, 0.05 * 0.151.^iter, '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'DisplayName', sprintf('Ref=$0.151$'));
plot(iter, 0.05 * 0.0062.^iter, '--', 'LineWidth', 2, 'Color', [0.6 0.4 0.4], 'DisplayName', sprintf('Ref=$0.0062$'));
ylim([1e-12,2]);
set(gca, 'YScale', 'log');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
xlabel('Iteration $k$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$e$', 'Interpreter', 'latex', 'FontSize', 20);
legend('Interpreter', 'latex', 'FontSize', 16, 'Location', 'northeast');
box on;
grid off;
hold off;

set(fig1, 'PaperPositionMode', 'auto');
set(fig2, 'PaperPositionMode', 'auto');
set(fig3, 'PaperPositionMode', 'auto');
set(fig4, 'PaperPositionMode', 'auto');
