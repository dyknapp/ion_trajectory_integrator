%% Get the mathieu parameters using scipy

% You need to add the python path with an environment that has scipy
% installed.  Use python 3.10
pyenv('Version', 'C:\Users\dk\miniconda3\envs\gp\python.exe', 'ExecutionMode', 'OutOfProcess');
pyrun('import scipy');
pyrun('import numpy as np');

%% Plot them

qs  = pyrun('qs = np.linspace(0, 1, 1000)', 'qs');
a0s = double(pyrun('a0s = scipy.special.mathieu_a(0.0, qs)', 'a0s'));
a1s = double(pyrun('a1s = scipy.special.mathieu_a(1.0, qs)', 'a1s'));
a2s = double(pyrun('a2s = scipy.special.mathieu_a(2.0, qs)', 'a2s'));
b1s = double(pyrun('b1s = scipy.special.mathieu_b(1.0, qs)', 'b1s'));
b2s = double(pyrun('b2s = scipy.special.mathieu_b(2.0, qs)', 'b2s'));
b3s = double(pyrun('b3s = scipy.special.mathieu_b(3.0, qs)', 'b3s'));

% figure;
% hold on
% plot(qs, a0s, '-k');
% plot(qs, b1s, '-r');
% 
% plot(qs, a1s, '-.k');
% plot(qs, b2s, '-.r');
% 
% plot(qs, a2s, '--k');
% plot(qs, b3s, '--r');
% hold off

%%

%%% Figure definition
fh = figure(...
    Name            = '',...
    WindowStyle     = 'docked',...
    WindowState     = 'normal',...
    Color           = 'w',...
    Units           = 'pixel',...
    ... Position        = [0 0 args.HeightResolution*args.AspectRatio args.HeightResolution],...
    PaperUnits      = 'centimeters',...
    PaperSize       = [9.0 9.0/((1+sqrt(5))/2)],...
    PaperPosition   = [0 0 9.0 9.0/((1+sqrt(5))/2)],...
    Resize          = 'on'...
    );


samples = 8192;
qs  = double(pyrun(sprintf('qs = np.linspace(0, 1,   %d)', samples), 'qs'));
as  = double(pyrun(sprintf('a_s = np.linspace(0, 0.5, %d)', samples), 'a_s'));
a0s = double(pyrun('a0s = scipy.special.mathieu_a(0.0, qs)', 'a0s'));
a1s = double(pyrun('a1s = scipy.special.mathieu_a(1.0, qs)', 'a1s'));
a2s = double(pyrun('a2s = scipy.special.mathieu_a(2.0, qs)', 'a2s'));
b1s = double(pyrun('b1s = scipy.special.mathieu_b(1.0, qs)', 'b1s'));
b2s = double(pyrun('b2s = scipy.special.mathieu_b(2.0, qs)', 'b2s'));
b3s = double(pyrun('b3s = scipy.special.mathieu_b(3.0, qs)', 'b3s'));

% figure;
% hold on
% plot(qs, -a0s, '-k');
% plot(qs, b1s, '-r');
% hold off

parameter_region = zeros(samples);
for q = 1:samples
    for a = 1:samples
        if -a0s(q) > as(a) && b1s(q) > as(a)
            parameter_region(q, a) = parameter_region(q, a) + 0.5;
        end
    end
end
for q = 1:samples
    for a = 1:samples
        if -a0s(q) > as(a) && b1s(q) > as(a)
            parameter_region(round(q/2), round(a/2)) = parameter_region(round(q/2), round(a/2)) + 1./8;
        end
    end
end

imagesc(qs, as, parameter_region')
set(gca,'YDir','normal');
xlim([0, 1])
ylim([0, 0.3])
map = [linspace(1,  30/256, samples); ...
       linspace(1,  47/256, samples); ...
       linspace(1, 151/256, samples); ]';
colormap(map)

hold on
plot(qs, [-a0s(qs <= 0.706) b1s(qs > 0.706)], '-k', LineWidth=2);
hold off

hold on
plot(qs/2, [-a0s(qs <= 0.706) b1s(qs > 0.706)]/2, '-k', LineWidth=2);
hold off

ax = gca;
ax.PlotBoxAspectRatio      = [((1+sqrt(5))/2) 1 1];
ax.Box                     = 'on';
ax.TickLabelInterpreter    = 'latex';
ax.FontName                = 'TimesNewRoman';
ax.FontSize                = 12;
ax.XMinorTick              = 'on';
ax.YMinorTick              = 'off';
ax.XGrid                   = 'off';
ax.YGrid                   = 'on';


