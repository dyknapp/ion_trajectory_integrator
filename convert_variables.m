%% Load data
clear
simion_path = "SIM001_data";
electrode_names = ["test_setup_assembly.pa15.patxt", ...
                   "test_setup_assembly.pa16.patxt", ...
                   "test_setup_assembly.pa17.patxt", ...
                   "test_setup_assembly.pa2.patxt", ...
                   "test_setup_assembly.pa4.patxt", ...
                   "test_setup_assembly.pa9.patxt", ...
                   "test_setup_assembly.pa13.patxt", ...
                   "test_setup_assembly.pa6.patxt", ...
                   "test_setup_assembly.pa14.patxt", ...
                   "test_setup_assembly.pa5.patxt", ...
                   "test_setup_assembly.pa10.patxt", ...
                   "test_setup_assembly.pa3.patxt", ...
                   "test_setup_assembly.pa11.patxt", ...
                   "test_setup_assembly.pa1.patxt", ...
                   "test_setup_assembly.pa12.patxt", ...
                  ];

start_line = 22;
loadanyways = false;
if ~exist('dimensions', 'var') || loadanyways
    addpath(simion_path)
    [potential_maps, is_electrode, dimensions] = ...
        readFile(electrode_names, start_line);
    potential_maps = potential_maps / 10000.0;
    is_electrode = logical(is_electrode);
end

%% Combine electrodes

einzel1 = potential_maps(1, :, :, :);
einzel2 = potential_maps(2, :, :, :);
einzel3 = potential_maps(3, :, :, :);

middle1 = potential_maps(10, :, :, :) + potential_maps(11, :, :, :);
middle2 = potential_maps(12, :, :, :) + potential_maps(13, :, :, :);
middle3 = potential_maps(14, :, :, :) + potential_maps(15, :, :, :);

endcap1 = potential_maps(5, :, :, :) + potential_maps(6, :, :, :);
endcap2 = potential_maps(4, :, :, :) + potential_maps(7, :, :, :);

rfelecs = potential_maps(8, :, :, :) + potential_maps(9, :, :, :);

%% Save the result

einzel1 = reshape(einzel1, dimensions);
einzel2 = reshape(einzel2, dimensions);
einzel3 = reshape(einzel3, dimensions);
endcap1 = reshape(endcap1, dimensions);
endcap2 = reshape(endcap2, dimensions);
rfelecs = reshape(rfelecs, dimensions);
middle1 = reshape(middle1, dimensions);
middle2 = reshape(middle2, dimensions);
middle3 = reshape(middle3, dimensions);

potential_maps = zeros([9 size(einzel3)]);
potential_maps(1, :, :, :) = einzel3;
potential_maps(2, :, :, :) = einzel2;
potential_maps(3, :, :, :) = einzel1;
potential_maps(4, :, :, :) = endcap1;
potential_maps(5, :, :, :) = endcap2;
potential_maps(6, :, :, :) = rfelecs;
potential_maps(7, :, :, :) = middle1;
potential_maps(8, :, :, :) = middle2;
potential_maps(9, :, :, :) = middle3;

save("decel_potential.mat", "potential_maps", "dimensions", "is_electrode", '-v7.3')

%% Check result

% to_plot = 1.0 * einzel1 + 2.0 * einzel2 + 3.0 * einzel3;
% to_plot = endcap1;
% to_plot = endcap2;
to_plot = rfelecs;
imagesc(squeeze(to_plot(round(dimensions(1)/2), :, :)))
axis equal




