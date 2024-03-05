%% PlotPaperGraph
% Plot a graph illustrating the trap frequencies for a single ion in the
% trap.

clear;
load('results.mat');

% Configure the figure
clf;
set(gcf, 'Units', 'centimeters');
pos = get(gcf, 'Position');
set(gcf, 'Position', [ pos(1) pos(2) 9 8 ], 'Color', 'w');

% Plot the results
pastelBlue = [112 146 190]/255;
plot(TrapQs, theoreticalRadialFreq*1e-6, '-', 'Color', 'k'); hold on;
plot(TrapQs, theoreticalZFreq*1e-6, '-', 'Color', pastelBlue);
plot(TrapQs, pseudoX*1e-6, 'o', 'Color', 'k', 'MarkerSize', 6, 'LineWidth', 1.5);
h1 = plot(TrapQs, pseudoZ*1e-6, 'o', 'Color', pastelBlue, 'MarkerSize', 6, 'LineWidth', 1.5);
grid on
set(gca, 'GridLineStyle', ':');

% and full rf
load('results_fullrf.mat');
plot(TrapQs, pseudoX*1e-6, 'x', 'Color', 'k', 'MarkerSize', 9, 'LineWidth', 1.5);
h2 = plot(TrapQs, pseudoZ*1e-6, 'x', 'Color', pastelBlue, 'MarkerSize', 9, 'LineWidth', 1.5);

xlabel('$q_y$', 'Interpreter', 'Latex', 'FontSize', 10);
set(get(gca, 'XAxis'), 'TickLabelInterpreter', 'Latex', 'FontSize', 10);
ylabel('$\omega_u/2\pi$ (MHz)', 'Interpreter', 'Latex', 'FontSize', 10);
set(get(gca, 'YAxis'), 'TickLabelInterpreter', 'Latex', 'FontSize', 10);

legend([h1 h2], {'pseudopotential', 'full rf'}, 'Interpreter', 'Latex', 'box', 'off', 'Location', 'northwest');

annotation('textbox', 'LineStyle', 'none', 'String', '$\omega_r$', 'Interpreter', 'Latex', 'FontSize', 11, 'Units', 'centimeters', 'Position', [ 4 4 1 1 ]);
annotation('textbox', 'LineStyle', 'none', 'String', '$\omega_z$', 'Interpreter', 'Latex', 'FontSize', 11, 'Units', 'centimeters', 'Position', [ 4 2 1 1 ], 'Color', pastelBlue);

set(gcf, 'Units', 'centimeters');
pos = get(gcf, 'Position');
w = pos(3); 
h = pos(4);
p = 0.01;
set(gcf,...
  'PaperUnits','centimeters',...
  'PaperPosition',[p*w p*h w h],...
  'PaperSize',[w*(1+2*p) h*(1+2*p)]);
set(gcf, 'Renderer', 'painters')
saveas(gcf, 'ver1.pdf');