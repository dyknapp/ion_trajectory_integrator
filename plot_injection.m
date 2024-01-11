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

linecolors = turbo(10);

% Variables to set beforehand
simion_path = "";
electrode_names = ["essentials_only_rf.patxt", ...
                   "essentials_only_dc1.patxt", ...
                   "essentials_only_dc2.patxt", ...
                   "essentials_only_dc3.patxt", ...
                   "essentials_only_dc4.patxt", ...
                   "essentials_only_dc5.patxt", ...
                  ];

start_line = 19;

loadanyways = false;
if ~exist('dimensions', 'var') || loadanyways
    addpath(simion_path)
    [potential_maps, is_electrode, dimensions] = ...
        readFile(electrode_names, start_line);
    potential_maps = potential_maps / 10000.0;
    is_electrode = logical(is_electrode);
end
potential_maps = reshape(potential_maps, [length(electrode_names) dimensions]);

d = 1.0;                % mm

start_time =  0.0;      % us
end_time   =  250.0;     % us

m = 2.0;                % amu (e.g. 2.0 would be roughly correct for H2+)
q = 1.0;                % atomic units

maxdist =  0.001;      % mm

%% Pre_injection

cx = round(dimensions(1)/2);
cy = round(dimensions(2)/2);
cz = round(dimensions(3)/2);
ys = d*(1:dimensions(2)) - cy;
zs = d*(1:dimensions(3)) - cz;
potential = tensorprod(potential_maps, [1.0, 0.0, 0.0, 0.0, 0.0, 5.0]', 1);
plot(zs, squeeze(potential(cx, cy, :)), 'color', linecolors(1,:), LineWidth=3)

hold on
potential = tensorprod(potential_maps, [1.0, 5.0, 0.0, 0.0, 0.0, 1.5]', 1);
plot(zs, squeeze(potential(cx, cy, :)), '-', 'color', linecolors(8,:), LineWidth=3)
hold off

hold on
potential = tensorprod(potential_maps, [1.0, 5.0, 0.0, 0.0, 0.0, 5.0]', 1);
plot(zs, squeeze(potential(cx, cy, :)), '--', 'color', linecolors(3,:), LineWidth=4)
hold off

for x = (-30:12:30)
    xline(x, '--k')
end
xlabel("Axial Displacement from Trap Center ($\mathrm{mm}$)", Interpreter="latex")
ylabel("DC Potential + Pseudopotential ($\mathrm{eV}$)", Interpreter="latex")

%%

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