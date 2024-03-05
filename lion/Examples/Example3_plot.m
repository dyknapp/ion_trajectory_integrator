%% Rigid Bodies:
% Plot a graph of the rigid body motion for the paper.

keyframeColor = [ 1 1 1 ] * 0.7;

% Calculate angle of rod
rodEnds = [rodAtoms(1).ID rodAtoms(end).ID];
rX = diff(x(rodEnds, :),1);
rY = diff(y(rodEnds, :),1);
rZ = diff(z(rodEnds, :),1);
rR = (rX.^2 + rY.^2).^0.5;
theta = rad2deg(atan(rZ./rR));
phi = wrapToPi(atan2(rY, rX));


% Setup figure size for paper.
clf
set(gcf, 'Color', 'w', 'Units', 'centimeters');
pos = get(gcf, 'Position');
set(gcf, 'Position', [ pos(1) pos(2) 9 9 ]);

py = 0.7;
% Create axes layout
ax_a = axes('Units', 'centimeters', 'Position', [ 0.35 0.6 2.6 3*2.6 ]); box on;
ax_theta = axes('Units', 'centimeters', 'Position', [ 4.1 6.4 4.7 2.0 ]); box on;
ax_phi = axes('Units', 'centimeters', 'Position', [ 4.1 3.9 4.7 2.0 ]); box on;

ax_b = axes('Units', 'centimeters', 'Position', [ 4.1 py 2.0 2.0 ]); box on;
ax_c = axes('Units', 'centimeters', 'Position', [ 6.7 py 2.0 2.0 ]); box on;
% ax_d = axes('Units', 'centimeters', 'Position', [ 3.4 10+2.4+py 2.4 2.4 ]); box on;
% ax_e = axes('Units', 'centimeters', 'Position', [ 6.4 10+2.4+py 2.4 2.4 ]); box on;
annotations_ax = axes('Position', [ 0 0 1 1 ]);
fp = get(gcf, 'Position');
xlim(annotations_ax, [ 0 fp(3) ]); ylim(annotations_ax, [ 0 fp(4) ]); hold(annotations_ax, 'on'); set(annotations_ax, 'Visible', 'off');

barc = [ 0.5 0.5 0.5 ];

% Subfigure a)
%  atoms at final positions.
axes(ax_a); frame = size(x,2);
h = depthPlot(x(:,frame),y(:,frame),z(:,frame-9), c, [ 20 70 ]*1.3, 'filled', [ 0.2 1 ], 'LineWidth', 0.5, 'MarkerEdgeColor', 'k');
view(0,0); view(-45,10);
% axis equal
set(gca, 'XTick', -525:25:500, 'XTickLabel', {});
set(gca, 'YTick', -525:25:500, 'YTickLabel', {});
set(gca, 'ZTick', -500:25:500, 'ZTickLabel', {});
xlim([-25 25]); ylim([-25 25]); zlim([-75 75]);
set(get(gca, 'XAxis'), 'TickDirection', 'in', 'TickLabelInterpreter', 'Latex', 'FontSize', 10);
set(get(gca, 'YAxis'), 'TickDirection', 'in', 'TickLabelInterpreter', 'Latex', 'FontSize', 10);
set(get(gca, 'ZAxis'), 'TickDirection', 'in', 'TickLabelInterpreter', 'Latex', 'FontSize', 10);
set(gca, 'GridLineStyle', '-');
box on;
l = xlabel(ax_b, 'x', 'Interpreter', 'Latex', 'Fontsize', 10); set(l, 'Units', 'centimeters', 'Position', [ 0.25 0.3 1 ]);
l = ylabel(ax_b, 'y', 'Interpreter', 'Latex', 'Fontsize', 10); set(l, 'Units', 'centimeters', 'Position', [ 1.6 0.3 1 ]);
hold on;
plot3([ 25 25; 25 25; 25 25; ]', -[ 0 0; 0 25; 25 25 ]', [ -69 -75; -72 -72; -69 -75 ]', 'Color', barc, 'LineWidth', 1.25);
tb = annotation('textbox', 'Interpreter', 'Latex', 'String', '$25 \mu$m', 'LineStyle', 'none', 'Units', 'centimeters', 'Position', [ 1.85 1.02 3 1 ], 'FontSize', 10);
set(get(tb, 'Text'), 'Rotation', -29)
h = xlabel(ax_a, 'x', 'Interpreter', 'Latex', 'FontSize', 10);
set(h, 'Units', 'centimeters', 'Position', [ 2 0.4 1 ]);
h = ylabel(ax_a, 'y', 'Interpreter', 'Latex', 'FontSize', 10);
set(h, 'Units', 'centimeters', 'Position', [ 0.6 0.4 1 ]);
h = zlabel(ax_a, 'z', 'Interpreter', 'Latex', 'FontSize', 10);
set(h, 'Units', 'centimeters', 'Position', [ 0 3.9 1 ]);

panelNum = 4; panelTimestepSpan = 2000;
frameNum = @(panel) floor(size(x,2) + (panel-1) / panelNum * panelTimestepSpan - panelTimestepSpan);

subfigs = [ ax_b ax_c ];
for i=1:2
% Subfigures (b)-(d)
%  Position of rod at the end of the simulation.
axes(subfigs(i)); frame = frameNum(i);
h = depthPlot(x(~rodMask,frame),y(~rodMask,frame),z(~rodMask,frame), c(~rodMask,:), [ 20 70 ], 'filled', [ 0.2 1 ], 'LineWidth', 0.5, 'MarkerEdgeColor', 'k', 'DistanceLimits', [ -75 75 ]); hold on;
set(h, 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0);
h = depthPlot(x(rodMask,frame),y(rodMask,frame),z(rodMask,frame), c(rodMask,:), [ 20 70 ], 'filled', [ 0.2 1 ], 'LineWidth', 0.5, 'MarkerEdgeColor', 'k', 'DistanceLimits', [ -75 75 ]); hold on;
view(0,90);
set(gca, 'XTick', -525:25:500, 'XTickLabel', {});
set(gca, 'YTick', -525:25:500, 'YTickLabel', {});
xlim([-25 25]); ylim([-25 25]);
set(get(gca, 'XAxis'), 'TickDirection', 'in', 'TickLabelInterpreter', 'Latex', 'FontSize', 10);
set(get(gca, 'YAxis'), 'TickDirection', 'in', 'TickLabelInterpreter', 'Latex', 'FontSize', 10);
set(gca, 'GridLineStyle', '-');
box on;


% Create scale bar aligned to graph divisions. Must be above axis.
p = get(gca, 'Position'); xl = xlim(gca); yl = ylim(gca);
u2yp = @(u) ((u-yl(1))/diff(yl)) .* p(4) + p(2);
u2xp = @(u) ((u-xl(1))/diff(xl)) .* p(3) + p(1);
barc = [ 0.5 0.5 0.5 ];
plot(annotations_ax, u2xp([ 0  0  ]), u2yp([ -22.5 -27.5 ]), 'Color', barc, 'LineWidth', 1.25);
plot(annotations_ax, u2xp([ 0  25 ]), u2yp([ -25.0 -25.0 ]), 'Color', barc, 'LineWidth', 1.25);
plot(annotations_ax, u2xp([ 25 25 ]), u2yp([ -22.5 -27.5 ]), 'Color', barc, 'LineWidth', 1.25);

end

ylabel(ax_b, 'y', 'Interpreter', 'Latex', 'FontSize', 10);
xlabel(ax_b, 'x', 'Interpreter', 'Latex', 'FontSize', 10);
xlabel(ax_c, 'x', 'Interpreter', 'Latex', 'FontSize', 10);
tb = annotation('textbox', 'Interpreter', 'Latex', 'String', '$25 \mu$m', 'LineStyle', 'none', 'Units', 'centimeters', 'Position', [ 7.6 0.2 3 1 ]);

% Subfigure f:
%  Angle of rod over time
t = timestep * sim.TimeStep * 1e3;
% Add lines for each sub diagram
xlim(ax_theta, [ min(t) max(t) ]);
ylim(ax_theta, [ -pi/2 pi/2 ]);
hold(ax_theta,'on');
plot(ax_theta, frameNum(1) * [1 1] * sim.TimeStep * 1e3, ylim, '-', 'Color', keyframeColor);
plot(ax_theta, frameNum(2) * [1 1] * sim.TimeStep * 1e3, ylim, '-', 'Color', keyframeColor);
plot(ax_theta, t, deg2rad(theta), 'k-');
set(get(ax_theta,'XAxis'), 'TickLabelInterpreter', 'Latex', 'FontSize', 10);
set(get(ax_theta,'YAxis'), 'TickLabelInterpreter', 'Latex', 'FontSize', 10);
set(ax_theta, 'GridLineStyle', ':');
set(ax_theta, 'YTick', [-pi/2 0 pi/2]); 
set(ax_theta, 'YTickLabels', {'$-\frac{\pi}{2}$','0', '$\frac{\pi}{2}$'});
set(ax_theta, 'XTick', [ 0 0.02 0.04 0.06 ], 'XTickLabels', {''});
l = ylabel(ax_theta, '$\theta$ (rad)', 'Interpreter', 'Latex', 'FontSize', 10);
set(l, 'Units', 'centimeters', 'Position', [ -0.6 1 0 ]);
grid(ax_theta,'on');

% Subfigure phi:
%  Angle phi of rod over time
t = timestep * sim.TimeStep * 1e3;
xlim(ax_phi, [ min(t) max(t) ]);
ylim(ax_phi, [ -pi pi ]);
hold(ax_phi,'on');
plot(ax_phi, frameNum(1) * [1 1] * sim.TimeStep * 1e3, ylim, '-', 'Color', keyframeColor);
plot(ax_phi, frameNum(2) * [1 1] * sim.TimeStep * 1e3, ylim, '-', 'Color', keyframeColor);
plot(ax_phi, t, phi, 'k-');
set(get(ax_phi,'XAxis'), 'TickLabelInterpreter', 'Latex', 'FontSize', 10);
set(get(ax_phi,'YAxis'), 'TickLabelInterpreter', 'Latex', 'FontSize', 10);
xlabel(ax_phi, 'Time (ms)', 'Interpreter', 'Latex', 'FontSize', 10);
set(ax_phi, 'GridLineStyle', ':');
set(ax_phi, 'YTick', [ -pi 0 pi ]); 
set(ax_phi, 'YTickLabels', {'$-\pi$', '0', '$\pi$'});
l = ylabel(ax_phi, '$\phi$ (rad)', 'Interpreter', 'Latex', 'FontSize', 10);
set(l, 'Units', 'centimeters', 'Position', [ -0.6 1 0 ]);
grid(ax_phi,'on');

% Axes behind plot to connect time frames to the associated preview boxes
keyframeColor = [ 0.9 0.9 0.9 ];
bgAx = axes('units', 'normalized', 'position', [ 0 0 1 1 ]);
uistack(bgAx, 'bottom');
v = 0.4;
xlim(bgAx, [ 0 1 ]); ylim(bgAx, [ 0 1 ]); hold(bgAx, 'on');
plot(0.508 * [ 1 1 ], [ 0.5 0.8 ], '-', 'Color', keyframeColor, 'LineWidth', 3);
plot(0.56 * [ 1 1 ], [ 0.5 0.8 ], '-', 'Color', keyframeColor, 'LineWidth', 3);
plot( [ 0.87 0.56 ], [ 0.295 0.432 ], '-', 'Color', keyframeColor, 'LineWidth', 3);
plot( [ 0.57 0.508 ], [ 0.295 0.432 ], '-', 'Color', keyframeColor, 'LineWidth', 3);
set(bgAx, 'Visible', 'off')


annotation('textbox', 'String', '(a)', 'FontSize', 11, 'LineStyle', 'none', 'Interpreter', 'Latex', 'Units', 'centimeters', 'Position', [ -0.15 8 1 1 ])
annotation('textbox', 'String', '(b)', 'FontSize', 11, 'LineStyle', 'none', 'Interpreter', 'Latex', 'Units', 'centimeters', 'Position', [ 2.9 8 1 1 ])
annotation('textbox', 'String', '(c)', 'FontSize', 11, 'LineStyle', 'none', 'Interpreter', 'Latex', 'Units', 'centimeters', 'Position', [ 2.9 2.05 1 1 ])

uistack(annotations_ax, 'top');

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
saveas(gcf, 'example3.pdf')