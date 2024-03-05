classdef InputFileElement < matlab.mixin.Heterogeneous & handle
    %INPUTFILEELEMENT Represents a portion of the input file
    %   Input file elements are used to construct the input file from
    %   modular pieces.
    
    properties        
        % A delegate that accepts varargin of the form (ID, InputFileArgs)
        createInputFileText;
        
        % A cell of arguments to be based to createInputFileText
        InputFileArgs;
    end
    
    properties (GetAccess=protected)
        % unique string ID representing this input file element. May be
        % used to reference it by other elements.
        ID;
    end
    
    methods
        
        function text = getLines(inputFileElement)
        %GETLINES returns the lines of text describing the given input
        %file element in LAMMPS commands.
            ife = inputFileElement;
            if iscell(ife.InputFileArgs)
                text = ife.createInputFileText(getID(ife), ife.InputFileArgs{:});
            elseif ~isempty(ife.InputFileArgs)
                text = ife.createInputFileText(getID(ife), ife.InputFileArgs);
            else
                text = ife.createInputFileText(getID(ife));
            end
        end
        
        function ref = getID(inputFileElement)
        %GETID Gets the unique ID representing the input file element.
            ref = inputFileElement.ID;
        end            
    end
    
end

