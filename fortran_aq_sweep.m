loadanyways = false;
if ~exist('dimensions', 'var') || loadanyways
    load("test_potential.mat");
end
electrode_names = ["einzel3" "einzel2" "einzel1" "endcap1" "endcap2" "rfelecs"];

RF_frequency = 13.2e+6;
d = 0.125;
m = 2 * 1.00727647 + 0.00054858;
q = 1.0;
end_time = 10.0;
maxdist  = 1.0e-2;

time_steps_per_us = 100;
time_steps = round((end_time - 0.0) * time_steps_per_us);
step_times = linspace(0.0, end_time + ((end_time - 0.0) / 100), time_steps);
voltages = zeros([time_steps 2]);

rf_electrode = 1;
endcaps = 2;

voltages(:, rf_electrode) = 1.0 * ...
    cos(2 * pi * step_times * RF_frequency / 10.0^6);
voltages(:, endcaps) = ones([time_steps 1]);

xx1 =      d * double(dimensions(1) + 1) / 2.0;
yy1 =      d * double(dimensions(2) + 1) / 2.0;
zz1 =      d * 400.0;

vxx1 = normrnd(0.0, 0.05);
vyy1 = normrnd(0.0, 0.05);
vzz1 = 10.0;

reps = 1;
res = 64;
a_res = res;
o_res = res;

amp_scales = linspace(0, 1200.0, a_res);
off_scales = linspace(-150, 150.0, o_res);

rf_potential = zeros([2 dimensions]);
rf_potential(1, :, :, :) = potential_maps(6, :, :, :);
endcap_potential = zeros([2 dimensions]);
endcap_potential(2, :, :, :) = 10.0 * (potential_maps(4, :, :, :) + potential_maps(5, :, :, :));

lifetimess = {};
parpool(4)
parfor vzz1 = double(1:30)
    tic
    lifetimes = fly_aqs(amp_scales, a_res, off_scales, o_res, reps, ...
                              xx1, yy1, zz1, vxx1, vyy1, vzz1, ...
                              rf_potential, endcap_potential, ... 
                              voltages, step_times, ...
                              time_steps, dimensions, int32(is_electrode), ...
                              int32(2), m, q, d, maxdist, end_time);
    newcell = {vzz1, lifetimes};
    lifetimess = [lifetimess; newcell];
    disp(vzz1)
    toc
end

%%

for idx = 1:30
    imagesc(amp_scales, off_scales, flip(lifetimess{idx, 2}', 1))
    set(gca,'YDir','normal')
    xlabel('RF amplitude (V)')
    ylabel('RF voltage offset (V)')
    title(sprintf('Beam energy / endcap potential = %0.2f', double(idx) / 22))
    exportgraphics(gcf,'testAnimated.gif','Append',true);
end
