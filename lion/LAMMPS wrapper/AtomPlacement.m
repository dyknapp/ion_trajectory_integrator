classdef AtomPlacement < handle
    %ATOMPLACEMENT Defines placement of an atom into simulation domain.
    
    properties (GetAccess=public, SetAccess=private)
        
        %X X-position of atom, m
        X
        
        %Y Y-position of atom, m
        Y
        
        %Z Z-position of atom, m
        Z
        
        %TYPE Atom type
        Type
        
        %ID ID associated with this atom
        ID
        
    end
    
    methods
        
        function instance = AtomPlacement(sim, x, y, z, type)
            %ATOMPLACEMENT Create an AtomPlacement instance.
            
            ip = inputParser;
            ip.addRequired('sim', @(sim) isa(sim, 'LAMMPSSimulation'));
            ip.addRequired('x', @(x) numel(x) == 1);
            ip.addRequired('y', @(y) numel(y) == 1);
            ip.addRequired('z', @(z) numel(z) == 1);
            ip.addRequired('type', @(type) isa(type, 'AtomType'));
            ip.parse(sim, x, y, z, type);
            
            instance.X = x;
            instance.Y = y;
            instance.Z = z;
            instance.Type = type;
            instance.ID = sim.AddAtom(instance);
        end
        
    end
    
end

