function [ g ] = Group( sim, content )
%GROUP Groups entities together in LAMMPS. 
% The body should either be an array of atomic species or an array of
% atomic ids. The first is assumed is input is numeric, the second if input
% is a struct with field 'id' specifying the atomic species id. If a match
% is found to an existing group the original is returned, otherwise a new
% group is created.
% 
% SYNTAX: sim.Group( 1:10 )
%         sim.Group( [ species1, species2 ] )

% Deduce group style from input type
if isa(content, 'AtomPlacement')
    style = 'id';
    c = [content.ID];
elseif isnumeric(content)
    style = 'id';
    c = content;
elseif isa(content, 'AtomType')
    style = 'type';
    c = [content.ID];
else
    error('Invalid content! Expect either an array of atomic ids or atomic species.');
end

% Create group from information
ng = LAMMPSGroup('style', style, 'content', c);

found = false;
% First see if another group exists matching this one.
for i = 1:length(sim.Groups)
    og = sim.Groups(i);
   
    if compare(og, ng)
        ng = og;
        found = true;
        break;
    end
end

if ~found
    sim.Groups(end+1) = ng;
end

g = ng;

end

