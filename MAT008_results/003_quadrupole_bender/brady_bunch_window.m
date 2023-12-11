function brady_bunch_window(ax, initial_location, side_lines, initial_spread, ...
                            initial_angle, initial_energy, voltages, ...
                            potential_maps, dimensions, is_electrode, electrode_names, ...
                            m, q, d, maxdist, end_time)
    constants = physical_constants();

    xx = initial_location(1);
    yy = initial_location(2);
    
    offsets = double(-side_lines:side_lines) * initial_spread * 0.5 / double(side_lines);
    particles = length(offsets);
    
    xxs = zeros([particles 1]);
    yys = zeros([particles 1]);
    
    for idx = 1:particles
        xxs(idx) = xx + offsets(idx) * -sin(initial_angle);
        yys(idx) = yy + offsets(idx) *  cos(initial_angle);
    end
    
    v = sqrt(2 * initial_energy * constants("elementary charge") / (m * constants("atomic mass unit"))) * 1.0e-3;
    vxxs = v * cos(initial_angle) * ones([particles 1]);
    vyys = v * sin(initial_angle) * ones([particles 1]);
    
    ms = m * ones([particles 1]);
    qs = q * ones([particles 1]);
    
    voltagess = repmat(voltages, [1, particles])';
    
    [x_trajs, y_trajs, tss, itss, datass] = ray_optics_ensemble(xxs, yys, vxxs, vyys, potential_maps, voltagess, dimensions, length(electrode_names), ms, qs, d, maxdist, end_time);
    %% 
    
    imagesc(d*(0:dimensions(1)-1), d*(1:dimensions(2)-1), is_electrode, 'Parent', ax)
    xlabel("x (mm)"); ylabel("y (mm)"); set(gca,'YDir','normal'); colormap("gray");
    axis image
    hold on
    for idx = 1:particles
        if offsets(idx) == 0
            ax.plot(x_trajs(idx, 1:datass(idx)), y_trajs(idx, 1:datass(idx)), '-r');
        else
            ax.plot(x_trajs(idx, 1:datass(idx)), y_trajs(idx, 1:datass(idx)), '-g', 'LineWidth', 0.1);
        end
    end
    hold off
end