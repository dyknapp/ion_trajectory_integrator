%% Initialize
constants = physical_constants();
fragment_energy_eV = 0.6;
dissoc_speed = 1.0e-3 * sqrt(2 * fragment_energy_eV * constants("elementary charge") ...
                                / constants("proton mass"));
batches = 24;
particles_per_batch = 256;
displacement = 3.0;

constants = physical_constants();
simion_path = "SIM001_data/005_TOF_electrode_stack";
electrode_names = ["backingplate.patxt", ...
                   "electrode1.patxt", ...
                   "electrode2.patxt", ...
                   "electrode3.patxt", ...
                   "electrode4.patxt", ...
                   "electrode5.patxt", ...
                   "electrode6.patxt", ...
                   "detector.patxt", ...
                  ];
start_line = 19;
loadanyways = true;
if ~exist('dimensions', 'var') || loadanyways
    addpath(simion_path)
    [potential_maps, is_electrode, dimensions] = ...
        readFile(electrode_names, start_line);
    % potential_maps = potential_maps / 10000.0;
    is_electrode = logical(is_electrode);
end

%% Optimize
% After a few hours of optimization:
% result =
% 
%  12.5829    0.1043    0.7949    1.0000    0.4009    0.1253    0.7664

guess = double([15, 10, 1, 1, 1, 3, 3]);

options = optimoptions( 'surrogateopt', ...
                        'Display','iter', ...
                        'PlotFcn', 'surrogateoptplot', ...
                        'MaxTime', Inf, ...
                        'ObjectiveLimit', 0.1, ...
                        'UseParallel', false);
options.MaxFunctionEvaluations = 1e6;
func = @(x)optimization_target( x, batches, particles_per_batch, ...
                                dissoc_speed, d, displacement, ...
                                potential_maps, is_electrode, dimensions, electrode_names);
result = surrogateopt(func, [0 10 -1000 -1000 -1000 -1000 -1000], [1 5000 0 0 0 0 0], options);

function score = optimization_target(x, batches, particles_per_batch, ...
                                     dissoc_speed, d, displacement, ...
                                     potential_maps, is_electrode, dimensions, electrode_names)
    n_particles = batches * particles_per_batch;
    xx1 = d * double(dimensions(1) + 1) / 2.0;  % mm
    yy1 = d * double(dimensions(2) + 1) * 0.91; % mm
    zz1 = d * double(dimensions(3) + 1) / 2.0;  % mm
    xxs = normrnd(xx1, displacement, [n_particles 1]);
    yys = normrnd(yy1, displacement, [n_particles 1]);
    zzs = normrnd(zz1, displacement, [n_particles 1]);
    score1 = ...
        score_stack_ours(x(1), ...
                    x(2), ...
                    x(3), ...
                    x(4), ...
                    x(5), ...
                    x(6), ...
                    batches, particles_per_batch, ...
                    dissoc_speed, ...
                    xxs, yys, zzs, ...
                    potential_maps, is_electrode, dimensions, electrode_names);
    score2 = ...
        score_stack_ours(x(2), ...
                    x(3), ...
                    x(4), ...
                    x(5), ...
                    x(6), ...
                    x(7), ...
                    batches, particles_per_batch, ...
                    dissoc_speed / 2., ...
                    xxs, yys, zzs, ...
                    potential_maps, is_electrode, dimensions, electrode_names);
    score4 = ...
        score_stack_ours(x(2), ...
                    x(3), ...
                    x(4), ...
                    x(5), ...
                    x(6), ...
                    x(7), ...
                    batches, particles_per_batch, ...
                    dissoc_speed / 4., ...
                    xxs, yys, zzs, ...
                    potential_maps, is_electrode, dimensions, electrode_names);
    score = score1 + score2 + score4;
    fprintf("SCORE: %14.8g  {R: %9.3g, W1: %9.3g, N1: %9.3g, W2: %9.3g, N2: %9.3g, W3: %9.3g}\n", ...
                score, x(2), x(3), x(4), x(5), x(6), x(7))
end