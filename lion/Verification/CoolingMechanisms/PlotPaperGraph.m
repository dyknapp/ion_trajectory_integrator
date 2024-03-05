%% PlotPaperGraph
% Plot a graph illustrating the dynamics of a Langevin bath.

clear;
load('results.mat');

% Configure the figure
clf;
set(gcf, 'Units', 'centimeters');
pos = get(gcf, 'Position');
set(gcf, 'Position', [ pos(1) pos(2) 9 8 ], 'Color', 'w');

% Get colors to use for plot. We interleave light and dark color for
% clarity.
pastelBlue = [112 146 190]/255;
black = [ 0 0 0 ];
colors = interp1( [ 0 1 ], [black; pastelBlue], linspace(0,1,size(T,2)));
colors = [ (colors(1:2:end,:)); colors(1:2:end,:); ];

% Plot temperatures and theory
for i=1:size(T,2)
    plot(t*1e6, T(:,i)*1e6, '-', 'Color', colors(i,:)); hold on;
    theoryT = (initialT + (bathT(i) - initialT)*(1-exp(1).^(-2*t/timeConstant(i))));
    plot(t*1e6, theoryT*1e6, ':', 'Color', colors(i,:));
end

grid on
set(gca, 'GridLineStyle', ':');

xlabel('$t$ ($\mu$s)', 'Interpreter', 'Latex', 'FontSize', 10);
set(get(gca, 'XAxis'), 'FontSize', 10, 'TickLabelInterpreter', 'Latex');
xlim([ 0 50 ]);

ylabel('$T_j$ ($\mu$K)', 'Interpreter', 'Latex', 'FontSize', 10);
set(get(gca, 'YAxis'), 'FontSize', 10, 'TickLabelInterpreter', 'Latex');
% ylim([ -3.5 3.5 ]);

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
saveas(gcf, 'ver4.pdf');