classdef LAMMPSGroup < handle
    %LAMMPSGROUP A group of atoms in LAMMPS
    
    properties
        %TYPE - the type of group.
        Style;
        
        %CONTENT - atom ids/styles that comprise the group
        Content;
        
        %ID - unique identifier of this group
        ID;
        
        %RIGID - flag specifying whether the group should be evolved as a
        %rigid body. Note that LION will not check whether atom assignments
        %to rigid groups overlap. It is assumed that a rigid group contains
        %atoms that are not in any other rigid group.
        Rigid;
    end
    
    methods
        
        function self = LAMMPSGroup(varargin)
            %LAMMPSGroup Creates a lammps group.
            % style: id or type
            % content: atom ids/types comprising group as a single vector.
            
            p = inputParser;
            
            validStyles = {'id', 'type'};
            p.addParameter('style','id', @(x) any(validatestring(x,validStyles)));
            p.addParameter('content',[], @(x) isnumeric(x));
            
            p.parse(varargin{:})
            
            self.Style = p.Results.style;
            self.Content = p.Results.content;
            
            self.ID = ['grp_' idString(5)];
            
            self.Rigid = 0;
            
        end
        
        function out = compare(grp1, grp2)
            %COMPARE compares two groups and returns true if grp1 = grp2
            % SYNTAX: compare(grp1, grp2)
            
            out = 0;
            
            % First, check style is the same
            if ~strcmpi(grp1.Style, grp2.Style)
                return;
            end
            
            % Second, check contents match
            if length(grp1.Content) ~= length(grp2.Content)
                return;
            end
            
            % comparison for floats
            if any(abs(sort(grp1.Content) - sort(grp2.Content)) > 0.001)
                return;
            end
            
            out = 1;
            
        end
        
        function ift = getInputFileText(self)
            %GETINPUTFILETEXT Gets input file text that defines the group.
            % SYNTAX: getInputFileText( group )
            
            content = sprintf('%d ', self.Content);
            %group groupName type content
            ift = sprintf('group %s %s %s\n', self.ID, self.Style, content);
        end
        
    end
    
end

