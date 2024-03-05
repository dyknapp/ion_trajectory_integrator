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

% Plot the predicted theoretical positions
pastelBlue = [112 146 190]/255;
for i=1:10
   h1 = plot(i, getTheoreticalPosition(i,1), '.k', 'MarkerSize', 10); hold on; 
end

for i=1:length(results)
    result = results(i);
    h2 = plot(result.N, result.positions/lengthScale, 'ok', 'Color', pastelBlue, 'LineWidth', 1.0, 'MarkerSize', 6);
end

grid on
set(gca, 'GridLineStyle', ':');

xlabel('$N$', 'Interpreter', 'Latex', 'FontSize', 10);
set(get(gca, 'XAxis'), 'FontSize', 10, 'TickLabelInterpreter', 'Latex');
set(gca, 'XTick', 1:10);
xlim([ 0.5 10.5 ]);

ylabel('$z/l$', 'Interpreter', 'Latex', 'FontSize', 10);
set(get(gca, 'YAxis'), 'FontSize', 10, 'TickLabelInterpreter', 'Latex');
set(gca, 'YTick', -3:1:3);
ylim([ -3.5 3.5 ]);

% Plot theory lines for odd and even numbers
% N = [results.N];
Nsmooth = linspace(0.5, 11, 500);
th = 2.018 ./ (Nsmooth).^0.559;
hold on;
plot(Nsmooth, th/2, 'k-')
plot(Nsmooth, -th/2, 'k-')

Nsmooth = linspace(0.5, 11, 500);
th = 2.018 ./ (Nsmooth).^0.559;
hold on;
plot(Nsmooth, th, 'k:')
plot(Nsmooth, 0*th, 'k:')
plot(Nsmooth, -th, 'k:')

legend([h1(1) h2(1)], {'theory', 'simulated'}, 'Interpreter', 'Latex', 'FontSize', 10, 'box', 'off', 'Units', 'centimeters', 'Position', [ 0.3 1.3 4 0.8 ]);

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
saveas(gcf, 'ver2.pdf');