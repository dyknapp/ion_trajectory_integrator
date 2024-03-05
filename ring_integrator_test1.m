
delta_r = 0.01;

z_separation = 1.0;
disc_r       = 10.0;

%% List of rings

ring_rs = [];
ring_zs = [];
ring_Vs = [];

% First disc
disc_voltage = 0.5;
disc_z       = z_separation / 2.;
for r = delta_r:delta_r:1
    ring_rs[length(ring_rs) + 1] = r;
    ring_zs[length(ring_zs) + 1] = disc_z;
    ring_Vs[length(ring_Vs) + 1] = disc_voltage;
end

% Second disc: only difference is negative signs for z and V
disc_voltage = -0.5;
disc_z       = -z_separation / 2.;
for r = delta_r:delta_r:1
    ring_rs[length(ring_rs) + 1] = r;
    ring_zs[length(ring_zs) + 1] = disc_z;
    ring_Vs[length(ring_Vs) + 1] = disc_voltage;
end

%% Construct the MoM matrix




