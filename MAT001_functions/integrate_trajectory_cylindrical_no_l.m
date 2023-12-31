function  output = integrate_trajectory_cylindrical_no_l(rr, zz, vrr, vzz, vOO, ...
                                        potential_maps, voltages, step_times, ...
                                        dimensions, is_electrode, ...
                                        m, q, d, maxdist, maxt)
% velocity verlet function with variable step size. only returns values of
% the end of the trajectory
%
% dknapp, 25.10.2023
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

    r = rr * 1.0e-3;
    z = zz * 1.0e-3;                % mm -> m
    v_r = vrr * 1.0e+3;
    v_z = vzz * 1.0e+3;             % mm/us -> m/s
    v_O = vOO;                      % already rad/s 
    step_times = step_times  * 1.0e-6;
    t = step_times(1);     % us -> s
    d = d * 1e-3;                   % mm -> m
    maxdist = maxdist * 1.0e-3;
    maxt = maxt * 1.0e-6;           % us -> s

    r_list = [r];
    z_list = [z];
    vr_list = [v_r];
    vz_list = [v_z];
    t_list = [t];
    
    % For cylindrical, we need the absolute m as well.
    m = (1.660539067e-27) * m;
    % Adjusting q / m units based on NIST numbers
    cmr = ((1.602176634e-19) * q) / m; % Coulomb/kg

    % reshape the potential maps in advance for speedy processing
    potential_maps_size = size(potential_maps);
    potential_maps = ...
        reshape(potential_maps, [potential_maps_size(1), (dimensions + int32([2, 0]))]);

    % Add the correction for the centripetal acceleration
    L = m * r^2 * v_O;

    [E_r, E_z] = field_at_cylindrical(r, z, d, potential_maps, ...
        interpolate_voltages(voltages, step_times(1), step_times));
    er_list = [E_r];
    ez_list = [E_z];

    % Also: it's a waste of time to calculate acceleration twice!!
    a_r = E_r * cmr;
    a_z = E_z * cmr;

    while t < maxt

        a = sqrt(a_r^2 + v_r^2);
        v = sqrt(v_r^2 + v_z^2);
      
        a = a + 1.0e-15;
        v = v + 1.0e-15;
        t_v = maxdist / v;
        t_a = sqrt(2 * maxdist / a);
        t_step = t_v * t_a / (t_v + t_a);

        t = t + t_step;
        
        % calculate coordinates for next timestep
        %  (m) + (s)*(m/s)    (s^2) * (kg m/s^2 C) * (C/kg)
        r = r + t_step * v_r + t_step^2 * a_r / 2;
        z = z + t_step * v_z + t_step^2 * a_z / 2;
        
        % stop loop if particle goes out of bounds
        % we need extra room so that the field calculation works
        if (z < 3 * d || ...
            r > double(dimensions(1) - 3) * d || z > double(dimensions(2) - 3) * d)
            break
        end

        if is_electrode(round(abs(r) / d) + 3, round(z/d))
            break
        end
        
        [E_r_new, E_z_new] = ...
            field_at_cylindrical(abs(r), z, d, potential_maps, ...
                     interpolate_voltages(voltages, t, step_times)); 

        % calculate velocities vor next timestep
        a_r_new = E_r_new * cmr;
        a_z_new = E_z_new * cmr;
        v_r = v_r + t_step * (a_r_new + a_r) / 2;
        v_z = v_z + t_step * (a_z_new + a_z) / 2;
            
        a_r = a_r_new;
        a_z = a_z_new;
        

        r_list(end + 1)  = r;
        z_list(end + 1)  = z;
        vr_list(end + 1) = v_r;
        vz_list(end + 1) = v_z;
        er_list(end + 1) = E_r;
        ez_list(end + 1) = E_z;
        t_list(end + 1)  = t;
    end
    
    %output coordinates
    output.r        = r_list * 1.0e+3;
    output.z        = z_list * 1.0e+3;
    output.vr       = vr_list * 1.0e-3;
    output.vz       = vz_list * 1.0e-3; % m / s -> mm/us
    output.er       = er_list * 1.0e-3; 
    output.ez       = ez_list * 1.0e-3; % V/mm -> V/m
    output.t        = t_list * 1.0e+6;  % s -> us
end