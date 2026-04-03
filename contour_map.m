% The plots show contour lines to indicate the peak shape and landscape.
% Fig 4.4

clear; close all; clc;
re = linspace(0.001, 10, 1000); 
im = linspace(-5, 5, 1000);
[Re, Im] = meshgrid(re, im);
S = Re + 1i * Im; 


gamma1 = zeros(size(S)); 
gamma2 = zeros(size(S)); 


denom = 1 + (2/3) * S;
R1 = (-1/3) ./ denom;
R2 = (4/3) ./ denom;
Nc = 1000;



for i = 1:numel(S)
    a = 1;
    b = -R2(i);
    c = -R1(i);
    disc = b^2 - 4*a*c;
    rho1 = (b + sqrt(disc)) / (2*a);
    rho2 = (b - sqrt(disc)) / (2*a);
    exp_neg_s = exp(-S(i));
    exp_neg_2s = exp_neg_s .^ 2;
    num = abs(exp_neg_2s - R2(i) * exp_neg_s - R1(i));
    % den = abs( (1 - abs(rho1)) * (1 - abs(rho2)) );

    a = rho1;
    b = rho2;
    diff_ab = abs(a - b);
        if diff_ab == 0
            den(i) = 0; 
        else
            sum_k = 0;
            for k = 0:Nc
                sum_k = sum_k + abs(a^(k+1) - b^(k+1)) / diff_ab;
            end
            den = sum_k;
        end
if den < 1e-10
        gamma1(i) = NaN;
else
        gamma1(i) = num * den;
end
end


denom = 1 + 0.56380 * S;
R1 = (0.02178 - 0.00047 * S) ./ denom;
R2 = (0.97822 - 0.46300 * S) ./ denom;



for i = 1:numel(S)
    a = 1;
    b = -R2(i);
    c = -R1(i);
    disc = b^2 - 4*a*c;
    rho1 = (b + sqrt(disc)) / (2*a);
    rho2 = (b - sqrt(disc)) / (2*a);
    exp_neg_s = exp(-S(i));
    exp_neg_2s = exp_neg_s .^ 2;
    num = abs(exp_neg_2s - R2(i) * exp_neg_s - R1(i));




if den < 1e-10
        gamma2(i) = NaN;
else
        gamma2(i) = num * den;
end
end

figure('Name', 'Contour Maps of \gamma^* for Re(s) > 0', 'NumberTitle', 'off', ...
'Position', [100, 100, 1200, 600], 'Color', 'w');



subplot(1,2,1);

levels = logspace(-3, 0, 21); 

[C, h] = contourf(Re, Im, gamma1, levels, 'LineWidth', 1.2);
clabel(C, h, 'FontSize', 12, 'FontName', 'Times New Roman', 'LabelSpacing', 200);
colormap(flipud(hot)); 
hold on;
plot(re, zeros(size(re)), 'k--', 'LineWidth', 0.5); 
plot(zeros(size(im)), im, 'k--', 'LineWidth', 0.5); 
xlim([min(re), max(re)]);
ylim([min(im), max(im)]);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman', ...
'TickLabelInterpreter', 'latex', 'Box', 'on', ...
'LineWidth', 1.2);
xlabel('$\Re(s)$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\Im(s)$', 'Interpreter', 'latex', 'FontSize', 18);

subplot(1,2,2);

levels = logspace(-3, 0, 10); 
levels = sort([levels, 0.0062]);
[C, h] = contourf(Re, Im, gamma2, levels, 'LineWidth', 1.2);
clabel(C, h, 'FontSize', 12, 'FontName', 'Times New Roman', 'LabelSpacing', 200);
colormap(flipud(hot));
hold on;
plot(re, zeros(size(re)), 'k--', 'LineWidth', 0.5); % Real axis
plot(zeros(size(im)), im, 'k--', 'LineWidth', 0.5); % Imaginary axis
xlim([min(re), max(re)]);
ylim([min(im), max(im)]);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman', ...
'TickLabelInterpreter', 'latex', 'Box', 'on', ...
'LineWidth', 1.2);
xlabel('$\Re(s)$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\Im(s)$', 'Interpreter', 'latex', 'FontSize', 18);

colorbar('FontSize', 14, 'LineWidth', 1.2, 'TickLabelInterpreter', 'latex');
