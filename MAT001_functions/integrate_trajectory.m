function  output = integrate_trajectory(xx, yy, zz, vxx, vyy, vzz, ...
                                        potential_maps, voltages, step_times, ...
                                        dimensions, is_electrode, ...
                                        m, q, d, maxdist, maxt)
% velocity verlet function with variable step size. only returns values of
% the end of the trajectory
%
% dknapp, 16.8.2023, based on TD_potential_TOF_VV3D.m
%
%
% INPUT:
%       xx, yy, zz          ... arrays of starting coordinates for N
%                               p   aricles in mm
%       vxx, vyy, vzz       ... arrays of starting velocities for N
%                                   paricles in mm/us
%       potential_maps      ... 3D array of potential values in V
%       step_times          ... Times for each step in parameter voltages
%       dimensions          ... dimensions of the potential array
%       z_detector          ... z coordinate in mm at which the detector
%                                   lies
%       m                   ... mass in amu
%       q                   ... charge in unitary charge units
%       d                   ... physical distance between potential gridpoints
%       maxdist             ... maximum distance a particle is allowed to
%                                   travel in mm
%       maxt                ... maximum time spent for the particle to
%                                   travel in microseconds
    

% OUTPUT:
%       x, y, z     ...     the trajectory
%       vx, vy, vz  ...     the velocities
%       t           ...     the time steps

    % f = waitbar(0,'1','Name','Integrating trajectory.', ...
    %     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

    % preallocate outputvariables

    x = xx * 1.0e-3;
    y = yy * 1.0e-3;
    z = zz * 1.0e-3;                % mm -> m
    v_x = vxx * 1.0e+3;
    v_y = vyy * 1.0e+3;
    v_z = vzz * 1.0e+3;             % mm/us -> m/s
    step_times = step_times  * 1.0e-6;
    t = step_times(1);     % us -> s
    d = d * 1e-3;                   % mm -> m
    maxdist = maxdist * 1.0e-3;
    maxt = maxt * 1.0e-6;           % us -> s

    x_list = [x];
    y_list = [y];
    z_list = [z];
    vx_list = [v_x];
    vy_list = [v_y];
    vz_list = [v_z];
    t_list = [t];

    % Adjusting q / m units based on NIST numbers
    cmr = ((1.602176634e-19) * q) / ((1.660539067e-27) * m); % Coulomb/kg

    % reshape the potential maps in advance for speedy processing
    potential_maps_size = size(potential_maps);
    potential_maps = ...
        reshape(potential_maps, [potential_maps_size(1), dimensions]);

    [E_x, E_y, E_z] = field_at(x, y, z, d, potential_maps, ...
        interpolate_voltages(voltages, step_times(1), step_times));
    ex_list = [E_x];
    ey_list = [E_y];
    ez_list = [E_z];

    while t < maxt
        a = sqrt(E_x^2 + E_y^2 + E_z^2) * abs(cmr);
        v = sqrt(v_x^2 + v_y^2 + v_z^2);
      
        a = a + 1.0e-15;
        v = v + 1.0e-15;
        t_v = maxdist / v;
        t_a = sqrt(2 * maxdist / a);
        t_step = t_v * t_a / (t_v + t_a);

        t = t + t_step;
        
        % calculate coordinates for next timestep
        %  (m) + (s)*(m/s)    (s^2) * (kg m/s^2 C) * (C/kg)
        x = x + t_step * v_x + t_step^2 * E_x * cmr / 2;
        y = y + t_step * v_y + t_step^2 * E_y * cmr / 2;
        z = z + t_step * v_z + t_step^2 * E_z * cmr / 2;
        
        % stop loop if particle goes out of bounds
        % we need extra room so that the field calculation works
        if (x < 3 * d || y < 3 * d || z < 3 * d || ...
            x > double(dimensions(1) - 3) * d || y > double(dimensions(2) - 3) * d || z > double(dimensions(3) - 3) * d)
            break
        end

        if is_electrode(round(x / d), round(y/d), round(z/d))
            break
        end
        
        [E_x_new, E_y_new, E_z_new] = ...
            field_at(x, y, z, d, potential_maps, ...
                     interpolate_voltages(voltages, t, step_times)); 

        % calculate velocities vor next timestep
        v_x = v_x + t_step * (E_x_new + E_x) * cmr / 2;
        v_y = v_y + t_step * (E_y_new + E_y) * cmr / 2;
        v_z = v_z + t_step * (E_z_new + E_z) * cmr / 2;
            
        E_x = E_x_new;
        E_y = E_y_new;
        E_z = E_z_new;
        
        x_list(end + 1) = x;
        y_list(end + 1) = y;
        z_list(end + 1) = z;
        vx_list(end + 1) = v_x;
        vy_list(end + 1) = v_y;
        vz_list(end + 1) = v_z;
        ex_list(end + 1) = E_x;
        ey_list(end + 1) = E_y;
        ez_list(end + 1) = E_z;
        t_list(end + 1) = t;
    end
    
    %output coordinates
    output.x        = x_list * 1.0e+3;
    output.y        = y_list * 1.0e+3;
    output.z        = z_list * 1.0e+3;  % m -> mm
    output.ex       = ex_list * 1.0e-3;
    output.ey       = ey_list * 1.0e-3;
    output.ez       = ez_list * 1.0e-3; % V/mm -> V/m
    output.vx       = vx_list * 1.0e-3;
    output.vy       = vy_list * 1.0e-3;
    output.vz       = vz_list * 1.0e-3; % m/s -> mm/us
    output.t        = t_list * 1.0e+6;  % s -> us
end