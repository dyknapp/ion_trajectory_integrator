res = 512;
dimensions = [res res];
potential_array = zeros(dimensions);
boundary_conditions = zeros(dimensions);

d = 100;

% Boundary conditions
for r = 1:dimensions(1)
    for z = 1:dimensions(2)
        % if    z >  round(dimensions(2) * 0.25) ...
        %    && z <  round(dimensions(2) * 0.75) ...
        %    && r == round(dimensions(1) * 0.27)
        %     boundary_conditions(r, z) = 10.0;
        % end
        % if    z >  round(dimensions(2) * 0.25) ...
        %    && z <  round(dimensions(2) * 0.75) ...
        %    && r == round(dimensions(1) * 0.25)
        %     boundary_conditions(r, z) = -1.0;
        % end
        if r == 3 * round(res / 10.0)
            if z > round(res / 4.0) && z < 3 * round(res / 4.0)
                boundary_conditions(r, z) = 10.0;
            end
        end
    end
end

imagesc(boundary_conditions)
set(gca,'YDir','normal')

for r = 2:dimensions(1)
    for z = 2:dimensions(2)
        if boundary_conditions(r, z) ~= 0.0
            potential_array(r, z) = boundary_conditions(r, z);
        end
    end
end

change = 0.0;
last_change = 0.0;
i = 0;
while (last_change - change) / (dimensions(1) * dimensions(2)) > 1.0e-10 ...
        || i <= 4 * max(dimensions)
    i = i + 1;
    last_change = change;
    change = 0.0;
    for r = 1:dimensions(1)
        potential_array(r, end) = potential_array(r, end - 1);
        potential_array(r, 1) = potential_array(r, 2);
    end
    for z = 1:dimensions(2)
        potential_array(end, z) = potential_array(end - 1, z);
    end
    for z = 2:dimensions(2) - 1
        for r = 2:dimensions(1) - 1
            % Non-boundaries
            if boundary_conditions(r, z) == 0.0
                prev = potential_array(r, z);
                potential_array(r, z) = (1 / 4.0) * ( ...
                     (d / (2.0 * r)) * (potential_array(r + 1, z) - potential_array(r + 1, z)) ...
                   + (potential_array(r + 1, z) + potential_array(r - 1, z)) ...
                   + (potential_array(r, z + 1) + potential_array(r, z - 1)));
                change = change + abs(potential_array(r, z) - prev);
            end
            potential_array(1, z) = (1 / 4.0) * (2 * potential_array(2, z) ...
                   + potential_array(1, z + 1) + potential_array(1, z - 1));
        end
    end
    fprintf("%d: %.3g, delta = %.3g\n", ...
        i, change / (dimensions(1) * dimensions(2)), ...
        (last_change - change) / (dimensions(1) * dimensions(2)))
end

% imcontour(potential_array(:, 2:end-1))
imagesc(potential_array(:, 2:end-1))
set(gca,'YDir','normal')