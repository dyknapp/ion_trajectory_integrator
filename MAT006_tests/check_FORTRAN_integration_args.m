function safe = check_FORTRAN_integration_args(xx, yy, zz, vxx, vyy, vzz, ...
                             potential_maps, voltages, step_times_in, ...
                             time_steps, dimensions, is_electrode, ...
                             n_electrodes, m, q, din, maxdist, maxt)
    safe = true;
    doubles_safe = false;

    try
        if ~isa(xx, 'double')
            fprintf("xx unsafe.\n");
            safe = false;
        elseif ~isa(yy, 'double')
            fprintf("yy unsafe.\n");
            safe = false;
        elseif ~isa(zz, 'double')
            fprintf("zz unsafe.\n");
            safe = false;
        elseif ~isa(vxx, 'double')
            fprintf("vxx unsafe.\n");
            safe = false;
        elseif ~isa(vyy, 'double')
            fprintf("vyy unsafe.\n");
            safe = false;
        elseif ~isa(vzz, 'double')
            fprintf("vzz unsafe.\n");
            safe = false;
        elseif ~isa(m, 'double')
            fprintf("m unsafe.\n");
            safe = false;
        elseif ~isa(q, 'double')
            fprintf("q unsafe.\n");
            safe = false;
        elseif ~isa(din, 'double')
            fprintf("din unsafe.\n");
            safe = false;
        elseif ~isa(maxdist, 'double')
            fprintf("maxdist unsafe.\n");
            safe = false;
        elseif ~isa(maxt, 'double')
            fprintf("maxt unsafe.\n");
            safe = false;
        end
        if safe
            doubles_safe = true;
        end
    %  real(c_double), dimension(time_steps, n_electrodes)
    % &      , intent(in) :: voltages
        if ~isequal(size(voltages), [time_steps n_electrodes])
            fprintf("voltages unsafe.\n");
            safe = false;
        end
    
    %  integer(c_int), dimension(3), intent(in) :: dimensions
        if ~isequal(size(dimensions), [1 3])
            fprintf("dimensions unsafe.\n");
            safe = false;
        end
    %  integer(c_int), intent(in),
    % &      dimension(dimensions(1), dimensions(2), dimensions(3))
    % &      :: is_electrode
        if ~isequal(size(is_electrode), dimensions)
            fprintf("is_electrode unsafe.\n");
            safe = false;
        end
    %  real(c_double), 
    % &      dimension(n_electrodes, dimensions(1), 
    % &      dimensions(2), dimensions(3)), 
    % &      intent(in) :: potential_maps
        if ~isequal(size(potential_maps), [n_electrodes dimensions])
            fprintf("potential_maps unsafe.\n");
            safe = false;
        end
    %  real(c_double), dimension(time_steps), intent(in) :: step_times_in
        if ~isequal(size(step_times_in), [1 time_steps])
            fprintf("step_times_in unsafe.\n");
            safe = false;
        end
    catch
        fprintf("Doubles safety: %d\n", doubles_safe)
        safe = false;
    end
end