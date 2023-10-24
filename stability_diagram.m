loadanyways = false;
if ~exist('dimensions', 'var') || loadanyways
    load("test_potential.mat");
end

res = 2048;
amps = linspace(0, 1200, res);
offs = linspace(-150, 150, res / 16);

success = zeros(res);

for a_idx = 1:res
    for q_idx = 1:res
        RF_amplitude = amps(q_idx);
        offset = offs(a_idx);
        m               = 2.0;
        q               = 1.0;
        RF_frequency    = 13.2e+6;
        
        
        end_time        = 5.0;
        maxdist         = 1.0e-3;
        
        vyy1            =  0.0;
        yy1             =  0.01;
        
        lifetime = try_real_trap_aq(RF_frequency, RF_amplitude, offset, 10.0, dimensions, is_electrode, potential_maps);
        fprintf("Ion lifetime in trap  was %.3g us (%d%% of the total simulation time)\n", ...
            lifetime, round(100.0 * lifetime / end_time))
        if lifetime >= end_time * 0.9
            success(a_idx, q_idx) = true;
        end
    end
end

%%
imagesc(amps, offs, flip(success(1:128, :), 1))
set(gca,'YDir','normal')
xlabel('RF amplitude (V)')
ylabel('RF voltage offset (V)')
title('Stability diagram for HD+, $f=13.2$MHz', 'Interpreter','latex')