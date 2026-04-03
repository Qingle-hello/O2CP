clc; clear; close all;

re = linspace(-20, 20, 1000);
im = linspace(-20, 20, 1000);
[Re, Im] = meshgrid(re, im);
S = Re + 1i * Im;


denom = 1 + (2/3) * S;
R1 = (-1/3) ./ denom;
R2 = (4/3) ./ denom;
a = 1;
b = -R2;
c = -R1;
disc = b.^2 - 4*a*c;
root1 = (b + sqrt(disc)) / (2*a);
root2 = (b - sqrt(disc)) / (2*a);
max_mod1 = max(cat(3, abs(root1), abs(root2)), [], 3);
stab1 = zeros(size(S));
stab1(max_mod1 <= 1 + 1e-10) = 1;

denom = 1 + 0.56380 * S;
R1 = (0.02178 - 0.00047 * S) ./ denom;
R2 = (0.97822 - 0.46300 * S) ./ denom;
a = 1;
b = -R2;
c = -R1;
disc = b.^2 - 4*a*c;
root1 = (b + sqrt(disc)) / (2*a);
root2 = (b - sqrt(disc)) / (2*a);
max_mod2 = max(cat(3, abs(root1), abs(root2)), [], 3);
stab2 = zeros(size(S));
stab2(max_mod2 <= 1 + 1e-10) = 1;

figure('Name', 'Stability Regions Comparison', 'NumberTitle', 'off', ...
'Position', [100, 100, 1200, 600], 'Color', 'w');

subplot(1,2,1);
contourf(Re, Im, stab1, [1 1], 'LineColor', 'none');
colormap([1 1 1; 0.8 0.9 1]);
hold on;
contour(Re, Im, stab1, [1 1], 'LineColor', 'k', 'LineWidth', 1.5);
contour(Re, Im, max_mod1, 0.1:0.1:1, 'LineColor', 'k', 'LineWidth', 1);
plot(re, zeros(size(re)), 'k--', 'LineWidth', 0.5);
plot(zeros(size(im)), im, 'k--', 'LineWidth', 0.5);
axis equal;
xlim([min(re), max(re)]);
ylim([min(im), max(im)]);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman', ...
'TickLabelInterpreter', 'latex', 'Box', 'on', ...
'LineWidth', 1.2);
xlabel('$\Re(s)$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\Im(s)$', 'Interpreter', 'latex', 'FontSize', 18);


subplot(1,2,2);
contourf(Re, Im, stab2, [1 1], 'LineColor', 'none');
colormap([1 1 1; 0.8 0.9 1]);
hold on;
contour(Re, Im, stab2, [1 1], 'LineColor', 'k', 'LineWidth', 1.5);
contour(Re, Im, max_mod2, 0.1:0.1:1, 'LineColor', 'k', 'LineWidth', 1);
plot(re, zeros(size(re)), 'k--', 'LineWidth', 0.5);
plot(zeros(size(im)), im, 'k--', 'LineWidth', 0.5);
axis equal;
xlim([min(re), max(re)]);
ylim([min(im), max(im)]);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman', ...
'TickLabelInterpreter', 'latex', 'Box', 'on', ...
'LineWidth', 1.2);
xlabel('$\Re(s)$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\Im(s)$', 'Interpreter', 'latex', 'FontSize', 18);

