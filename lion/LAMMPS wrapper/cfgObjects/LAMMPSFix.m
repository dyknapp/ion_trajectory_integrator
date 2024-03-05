classdef LAMMPSFix < InputFileElement
    %FIX All fixes
    
    properties
        %TIME apply a limiting timestep if time != 0;
        time
        
        %BASENAME prepend this string_ to fixIDs
        baseName
    end
    
    methods
        %Create a fix object, assigning it an unused ID.
        function obj = LAMMPSFix(baseName)
            if nargin < 1
                baseName = '';
            end
            obj.time = 0;
            obj.baseName = baseName;
        end
        
        function id = getID(obj)
           
            if ~isempty(obj.baseName)
                id = [obj.baseName '_' obj.ID];
            else
                id = obj.ID;
            end
            
        end
    end
end
