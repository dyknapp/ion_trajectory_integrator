classdef AtomType < handle
    %ATOMTYPE A type of atom in a simulation
    
    properties (SetAccess=private, GetAccess=public)
        
        %ID ID assigned to this atom type in the simulation
        ID;
        
        %CHARGE Charge of this atom type, in elementary charge units +e.
        Charge;
        
        %MASS Mass of this atom type, in atomic mass units amu.
        Mass
        
        %GROUP Group selection for this atom type.
        Group
        
    end
    
    methods
       
        function instance = AtomType(sim, id, charge, mass)
            %ATOMTYPE Creates a new atom type
            instance.ID = id;
            instance.Charge = charge;
            instance.Mass = mass;
            instance.Group = sim.Group(instance);
        end
        
    end
    
end

