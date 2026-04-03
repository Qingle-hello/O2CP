clc; clear; close all;
% Define iteration vector
iter_max = 10;
iter = 0:iter_max;


% Load data for nonsmooth, J=20
load('BDF2_solver_J_20.mat', 'err');
err_BDF2_J20_ns = err;
load('MOCP1_solver_J_20.mat', 'err');
err_MOCP1_J20_ns = err;
load('OCP_J_20.mat', 'err');
err_OCP_J20_ns = err;
load('SDIRK2_J_20.mat', 'err');
err_SDIRK2_J20_ns = err;
% Load data for nonsmooth, J=50
load('BDF2_solver_J_50.mat', 'err');
err_BDF2_J50_ns = err;
load('MOCP1_solver_J_50.mat', 'err');
err_MOCP1_J50_ns = err;
load('OCP_J_50.mat', 'err');
err_OCP_J50_ns = err;
load('SDIRK2_J_50.mat', 'err');
err_SDIRK2_J50_ns = err;


% First figure: Problem data (c), J=20
fig1 = figure('Position', [100, 100, 600, 400]);
hold on;
% Plot data lines with default colors, different markers, linewidth 3
plot(iter, err_BDF2_J20_ns, '-o', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'BDF2 ');
plot(iter, err_MOCP1_J20_ns, '-s', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'O2CP');
plot(iter, err_OCP_J20_ns, '-^', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'OCP');
plot(iter, err_SDIRK2_J20_ns, '-d', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'SDIRK2');
% Plot reference lines with dashed lines
plot(iter, 0.05 * 0.151.^iter, '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'DisplayName', sprintf('Ref=$0.151$'));
plot(iter, 0.05 * 0.0062.^iter, '--', 'LineWidth', 2, 'Color', [0.6 0.4 0.4], 'DisplayName', sprintf('Ref=$0.0062$'));
ylim([1e-12,2]);
% Set log scale for y-axis
set(gca, 'YScale', 'log');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
% Labels with LaTeX interpreter
xlabel('Iteration $k$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$e$', 'Interpreter', 'latex', 'FontSize', 20);
% title('Problem data (c), $J=20$', 'Interpreter', 'latex', 'FontSize', 20);
% Legend with LaTeX
legend('Interpreter', 'latex', 'FontSize', 16, 'Location', 'northeast');
% Box on, grid off, font settings
box on;
grid off;
hold off;
% Second figure: Problem data (c), J=50
fig2 = figure('Position', [700, 100, 600, 400]);
hold on;
plot(iter, err_BDF2_J50_ns, '-o', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'BDF2 ');
plot(iter, err_MOCP1_J50_ns, '-s', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'O2CP');
plot(iter, err_OCP_J50_ns, '-^', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'OCP');
plot(iter, err_SDIRK2_J50_ns, '-d', 'LineWidth', 3, 'MarkerSize', 12, 'DisplayName', 'SDIRK2');
% Plot reference lines with dashed lines
plot(iter, 0.05 * 0.151.^iter, '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'DisplayName', sprintf('Ref=$0.151$'));
plot(iter, 0.05 * 0.0062.^iter, '--', 'LineWidth', 2, 'Color', [0.6 0.4 0.4], 'DisplayName', sprintf('Ref=$0.0062$'));
ylim([1e-12,2]);
% Set log scale for y-axis
set(gca, 'YScale', 'log');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
xlabel('Iteration $k$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$e$', 'Interpreter', 'latex', 'FontSize', 20);
% title('Problem data (c), $J=50$', 'Interpreter', 'latex', 'FontSize', 20);
legend('Interpreter', 'latex', 'FontSize', 16, 'Location', 'northeast');
box on;
grid off;
hold off;
% Tighten layout
set(fig1, 'PaperPositionMode', 'auto');
set(fig2, 'PaperPositionMode', 'auto');
