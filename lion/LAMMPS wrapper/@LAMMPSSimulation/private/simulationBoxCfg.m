function [ strings ] = simulationBoxCfg( l, w, h, ionSpecies )
%SIMULATIONBOXCFG Generates cfg file text for simulation box.

if isempty(l) || size(l,1) == 0 || isempty(w) || size(w,1) == 0  || isempty(h) || size(h,1) == 0 
   error('Simulation box size undefined. Use experiment.SetSimulationDomain(l,w,h) to define the simulation region.'); 
end

strings = { '#Creating simulation box...' };
strings{end+1} = 'boundary m m m';
strings{end+1} = sprintf('region simulationDomain block %e %e %e %e %e %e units box', -l, l, -w, w, -h, h);
strings{end+1} = sprintf('create_box %d simulationDomain', ionSpecies); %<-- need to change this to include number of atom types in simulation!
strings{end+1} = '';
end

