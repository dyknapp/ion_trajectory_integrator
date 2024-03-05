classdef LAMMPSVariable < InputFileElement
    %LAMMPSVARIABLE Interface, represents some variable used for output or
    %further operation in lammps, eg average velocity, secular velocity,
    %etc.
    
    properties
        Frequency %function handle that returns the frequency in steps between when the variable is defined, eg for step 100, 200, 300... this value would be 100.
        Length
        Text      %used for atomProperty to specify which property to use.
    end
    
    properties (SetAccess = private)
        variableType
    end
    
    methods
        function obj = LAMMPSVariable()
            %LAMMPSVARIABLE Creates a new LAMMPSVariable object and assigns
            %it an unused id
            obj.ID =  idString(6);
            obj.Frequency = @() 1;
            obj.setVariableType('fix');
            obj.createInputFileText = @(~) { };
            obj.Length = 1;
        end
        
        function outputs = getOutputs(obj)
            %GETOUTPUTS Returns a cell containing string representations of
            %valid references to all outputs of this fix/variable in LAMMPS
            
            % Eg. returns { "f_a[1]", "f_a[2]" } for length-2 fix with ID
            % 'a', or { 'x' } for atomProperty x.
            
            switch obj.variableType
                case {'fix', 'var'}
                    base = [getVariablePrefix(obj.variableType) obj.ID];
                    outputs = cell(1,obj.Length);
                    if obj.Length > 1
                        for i = 1:obj.Length
                            outputs{i} = [base '[' num2str(i, '%d') ']' ];
                        end
                    else
                        outputs{1} = base;
                    end
                case 'atomProperty'
                    outputs = cell(1,1);
                    outputs{1} = obj.Text;
                otherwise
                    error('Invalid variable type.');
                    
            end
        end
        
        function setVariableType(obj, type)
            %SETVARIABLETYPE Assigns a type to variable. Throws an error for
            %invalid types.
            switch type
                case 'fix'
                case 'var'
                case 'atomProperty'
                    
                otherwise
                    error('Invalid variable type.');
                    
            end
            obj.variableType = type;
        end
    end
    
end

