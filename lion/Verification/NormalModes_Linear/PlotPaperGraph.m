%% PlotPaperGraph
% Plot a graph showing the normal mode spectra of the linear chain.

load('results.mat');

clf;
set(gcf, 'Color', 'w', 'Units', 'centimeters');
pos = get(gcf, 'Position');
set(gcf, 'Position', [ pos(1) pos(2) 9 8 ]);

index = 1;
xTick = [];
xTickStr = {};
for i=1:length(results)
    result = output(i);
    % plot a strip for each normal mode spectrum
    for j=1:result.N
        a = result.amps(:,j);
        
        % the amplitudes must be averaged, otherwise tiny features (such as
        % the resonances themselves) are lost to aliasing. Also saves on
        % image size.
        res = 10;
        a = movmean(a, res); a = a(1:res:end);
        f = movmean(result.f, res); f = f(1:res:end);
        
        
        span = [index-0.5 index+0.5 ]+i*0.5;
        x = repmat(span, length(a), 1);
        y = repmat(f'/(axialAngularFreq/(2*pi)), 1, 2);
        c = repmat(a, 1, 2);        
        surf(x, y, c, log(c/max(c(:)))); shading interp; hold on;
        
        % add lines top/bottom, otherwise axis is occluded
        ylim([ 0 10 ]);
%         plot3(span, [1 1]*min(ylim), [1 1], 'k-');
%         plot3(span, [1 1]*min(ylim), [1 1], 'k-');
        
        xTick(end+1) = index+i*0.5;
        xTickStr{end+1} = sprintf('%d', i);
        
        index = index+1;
        
        % plot theoretical points
        plot3(mean(span), result.freqs(j), 1, '.k');
    end
end

xlim([ 1 index-0.5+i*0.5 ]);

% Add dividers to separate different values of N
index = 0;
for i=1:length(results)
    plot3([ 1 1 ]*(index+0.5)+i*0.5, ylim, [1 1],'k-');
    plot3([ 1 1 ]*(index+0.5+results(i).N)+i*0.5, ylim, [1 1],'k-');
    index = index + results(i).N;
end

% Re-rule axes, which seems to be occluded by the surf plots
plot3(xlim, [1 1]*min(ylim), [1 1], 'k-');
plot3(xlim, [1 1]*max(ylim), [1 1], 'k-');
plot3([1 1]*min(xlim), ylim, [1 1], 'k-');
plot3([1 1]*max(xlim), ylim, [1 1], 'k-');

view([ 0 90 ]);

set(gca, 'XTick', [ results.N ].^2/2+[results.N]/2+0.5, 'XTickLabels', arrayfun(@(x) sprintf('%d', x), [ results.N ], 'UniformOutput', 0));
colormap(parulawhite);
grid off
caxis([ -14 0 ]);

xlabel('$N$', 'Interpreter', 'Latex', 'FontSize', 10);
set(get(gca, 'XAxis'), 'TickLabelInterpreter', 'Latex', 'FontSize', 10);

ylabel('$\omega_p / \omega_z$', 'Interpreter', 'Latex', 'FontSize', 10);
set(get(gca, 'YAxis'), 'TickLabelInterpreter', 'Latex', 'FontSize', 10, 'TickDirection', 'out');

print('ver3.png', '-r500', '-dpng');