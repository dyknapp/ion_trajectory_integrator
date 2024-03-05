%% Plot the results
% This plots the figure used in the paper.

[timestep, ~, x,y,z] = readDump('positions.txt');
x = x*1e6; y=y*1e6; z=z*1e6; % convert to um

% Setup figure size for paper.
clf
set(gcf, 'Color', 'w', 'Units', 'centimeters');
pos = get(gcf, 'Position');
set(gcf, 'Position', [ pos(1) pos(2) 9 7.5 ]);

% Subfigure a)
%  Illustration of trajectories along z axis of select atoms
pastelBlue = [112 146 190]/255;
ax = axes('Units', 'centimeters', 'Position', [ 1.1 3.95 7.5 3 ]);
plot(timestep(1:20:end), z(1:30,1:20:end)', 'Color', pastelBlue); hold on;
plot(timestep, z(1,:)', 'Color', [ 0.1 0.1 0.1 ], 'LineWidth', 1.0); hold off;
xlim([ 0 size(z, 2) ]); ylim([-150 150]); set(gcf, 'Color', 'w');
xlabel('timestep', 'Interpreter', 'Latex', 'FontSize', 10);
set(gca, 'XTick', [ 0, 5000, 10000, 15000 ], 'GridLineStyle', ':');
set(get(gca, 'YAxis'), 'TickDirection', 'in', 'TickLabelInterpreter', 'Latex', 'FontSize', 10);
set(get(gca, 'XAxis'), 'TickDirection', 'in', 'TickLabelInterpreter', 'Latex', 'FontSize', 10, 'Exponent', 3);
set(gca, 'YTick', -150:50:100);
grid on
lab = ylabel('z ($\mu$m)', 'Interpreter', 'Latex', 'FontSize', 10);
set(lab, 'Units', 'centimeters', 'Position', [ -0.7 1.5 0]);
get(lab, 'Position');


% Subfigure b)
%  Final trajectories of ions, resembling a Coulomb crystal. The crystal is
%  rendered in 3D.
ax = axes('Units', 'centimeters', 'Position', [ 1.1 0.5 2.2 2.2 ]);
depthPlot(x(:,end), z(:,end), y(:,end), pastelBlue, [ 20 70 ], 'filled', [ 0.2 1 ], 'LineWidth', 0.5, 'MarkerEdgeColor', 'k');
xlab = xlabel('x', 'Interpreter', 'Latex', 'FontSize', 10, 'Units', 'centimeters', 'Rotation', -0);
zlab = zlabel('y', 'Interpreter', 'Latex', 'FontSize', 10, 'Units', 'centimeters');
set(gca,'TickLength',[0.02 0.02], 'FontSize', 10);
set(get(gca, 'XAxis'), 'TickLabelInterpreter', 'Latex');
set(get(gca, 'YAxis'), 'TickLabelInterpreter', 'Latex');
set(get(gca, 'ZAxis'), 'TickLabelInterpreter', 'Latex');
grid on;  view(0,0); axis equal ;%axis vis3d
xlim([-50 50])
ylim([-150 150]);
zlim([-50 50]);
set(gca, 'GridLineStyle', '-');
box on
set(gca, 'XTick', -50:50:50, 'XTickLabel', {});
set(gca, 'YTick', -200:50:200, 'YTickLabel', {});
set(gca, 'ZTick', -50:50:50, 'ZTickLabel', {});

% Panel c)
ax = axes('Units', 'centimeters', 'Position', [ 4.2 0.5 4.4 2.2 ]);
depthPlot(x(:,end), z(:,end), y(:,end), pastelBlue, [ 20 70 ], 'filled', [ 0.2 1 ], 'LineWidth', 0.5, 'MarkerEdgeColor', 'k');
ylab = ylabel('z', 'Interpreter', 'Latex', 'FontSize', 10, 'Units', 'centimeters', 'Rotation', 0);
zlab = zlabel('y', 'Interpreter', 'Latex', 'FontSize', 10, 'Units', 'centimeters');
set(gca,'TickLength',[0.02 0.02], 'FontSize', 10);
set(get(gca, 'XAxis'), 'TickLabelInterpreter', 'Latex');
set(get(gca, 'YAxis'), 'TickLabelInterpreter', 'Latex');
set(get(gca, 'ZAxis'), 'TickLabelInterpreter', 'Latex');
grid on;  view(90,0); axis equal ;%axis vis3d
xlim([-50 50]);
ylim([-100 100]);
zlim([-50 50]);
set(get(gca, 'XAxis'), 'TickDirection', 'in', 'FontSize', 11)
set(get(gca, 'YAxis'), 'TickDirection', 'in')
set(get(gca, 'ZAxis'), 'TickDirection', 'in')
set(gca, 'GridLineStyle', '-');
box on
set(gca, 'XTick', -50:50:50, 'XTickLabel', {});
set(gca, 'YTick', -200:50:200, 'YTickLabel', {});
set(gca, 'ZTick', -50:50:50, 'ZTickLabel', {});
hold on


% Create scale bar aligned to graph divisions. Must be above axis.
p = get(ax, 'Position'); fp = get(gcf, 'Position'); xl = ylim(ax); yl = zlim(ax);
u2yp = @(u) ((u-yl(1))/diff(yl)) .* p(4) + p(2);
u2xp = @(u) ((u-xl(1))/diff(xl)) .* p(3) + p(1);
barc = [ 0.5 0.5 0.5 ];
ax = axes('Position', [ 0 0 1 1 ]);
plot(u2xp([ 50  50  ]), u2yp([ -45 -55 ]), 'Color', barc, 'LineWidth', 1.25); hold on;
plot(u2xp([ 50  100 ]), u2yp([ -50 -50 ]), 'Color', barc, 'LineWidth', 1.25); hold on;
plot(u2xp([ 100 100 ]), u2yp([ -45 -55 ]), 'Color', barc, 'LineWidth', 1.25); hold off;
xlim([ 0 fp(3) ]); ylim([ 0 fp(4) ]);
set(ax, 'Visible', 'off');
tb = annotation('textbox', 'Interpreter', 'Latex', 'String', '$50 \mu$m', 'LineStyle', 'none', 'Units', 'centimeters', 'Position', [ 7.45 -0.5 3 1 ]);


% Subfigure labels
annotation('textbox', 'String', '(a)', 'FontSize', 11, 'LineStyle', 'none', 'Position', [ -0.01 0.93 0.05 0.05 ], 'Interpreter', 'Latex')
annotation('textbox', 'String', '(b)', 'FontSize', 11, 'LineStyle', 'none', 'Position', [ -0.01 0.36 0.05 0.05 ], 'Interpreter', 'Latex')
annotation('textbox', 'String', '(c)', 'FontSize', 11, 'LineStyle', 'none', 'Position', [ 0.37 0.36 0.05 0.05 ], 'Interpreter', 'Latex')

% Render to file
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
saveas(gcf, 'example1.pdf')