classdef LAMMPSSimulation < handle & matlab.mixin.SetGet
    % LAMMPSSIMULATION A LAMMPS simulation that may be configured by defining atoms types, adding atoms and forces.
    
    properties
        
        %TIMESTEP Time between adjacent simulation steps in seconds
        TimeStep
        
        %CONFIGFILENAME Name of the generated LAMMPS input file.
        ConfigFileName
        
        %COULOMBCUTOFF Cut-off range for the Coulomb interaction
        CoulombCutoff
        
        %GPUACCEL Enable/Disable gpu acceleration through the gpu package.
        GPUAccel
        
        %NEIGHBORLIST Select style used for neighbour list
        NeighborList % (see http://lammps.sandia.gov/doc/neighbor.html)
        
        %NEIGHBORSKIN Size of the skin used for neighbor calculations
        NeighborSkin % (see http://lammps.sandia.gov/doc/neighbor.html)
        
        %VERBOSE Whether lammps terminal should be flushed to see progress.
        Verbose = 1;
        
    end
    
    properties (SetAccess=private)
        
        %ELEMENTS Time-Ordered list of simulation elements
        Elements
        
        %ATOMLIST list of function handles that will create atoms
        AtomList
        
        %SIMULATIONBOX Defines the simulation region
        SimulationBox
        
        %LIMITINGTIMESTEP Allows fixes which require a max timestep to
        %limit the timestep size.
        LimitingTimestep
        
        %ATOMTYPES List of different types of atom in the simulation
        AtomTypes
        
        %HASEXECUTED Has the simulation been run?
        HasExecuted
        
    end
    
    properties (Access=private)
        
        %GROUPS A list of groups defined in the simulation.
        Groups
        
    end
    
    methods
        
        function sim = LAMMPSSimulation
            % LAMMPSSimulation Create a simulation template for lammps
            % SYNTAX: LAMMPSSimulation()
            
            sim.Elements = InputFileElement.empty(1,0); %create empty array of fix objects.
            sim.AtomList = AtomPlacement.empty(1,0);
            sim.AtomTypes = AtomType.empty(1,0);
            sim.Groups = LAMMPSGroup.empty(1,0);
            sim.SimulationBox = struct('width', {}, 'height', {}, 'length', {});
            sim.TimeStep = 1;
            sim.LimitingTimestep = 1;
            sim.ConfigFileName = 'sim.lammps';
            sim.HasExecuted = false;
            sim.CoulombCutoff = 10;
            sim.GPUAccel = 0;
            sim.NeighborList = 'nsq';
            sim.NeighborSkin = 1;
        end
        
        SetSimulationDomain(obj, length, width, height)
        id = AddAtom(obj, atoms)
        [atomTypeStruct] = AddAtomType(obj, charge, mass);
        Execute(obj)
        WriteInputFile(obj)
        Unfix(obj, fix)
        Add(sim, obj)
        indices = GetSpeciesIndices(obj)
        Remove(sim, obj)
        g = Group(sim, content)
        
        function assertNotRun(sim)
            %ASSERTNOTRUN Asserts that the simulation has not yet run.
            if sim.HasExecuted
                error('Simulation has already executed. State cannot be modified.');
            end
        end
        
    end
end

