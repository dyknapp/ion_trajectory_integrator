res = 1000;
as = linspace(-0.25, 0.75, res);
qs = linspace(0, 1, res);

success = zeros(res);

parfor a_idx = 1:res
    for q_idx = 1:res
        aM = as(a_idx);
        qM = qs(q_idx);
        m               = 2.0;
        q               = 1.0;
        RF_frequency    = 10.0e+6;
        
        
        end_time        = 100.0;
        maxdist         = 1.0e-4;
        
        vyy1            =  0.0;
        yy1             =  0.01;
        
        result = try_mathieu_aq(aM, qM, m, q, RF_frequency, end_time, maxdist, yy1, vyy1);
        fprintf("Ion lifetime in trap  was %.3g us (%d%% of the total simulation time)\n", ...
            result.lifetime, round(100.0 * result.lifetime / end_time))
        if result.lifetime >= end_time - 0.001
            success(a_idx, q_idx) = true;
        elseif result.its == 2^20
            fprintf("Simulation was terminated early due to max timestep limit (2^20) being reached.\n")
            success(a_idx, q_idx) = true;
        end
    end
end