function [ runCmd ] = custom( command )
%CUSTOM Allows insertion of a custom script command into LAMMPS.
% 
% This command is intended as a way to allow a more advanced user to easily
% extend the wrapper's ability beyond the standard methods. Text entered
% using the command 'custom' will be inserted into the generated lammps
% input file.
%
% SYNTAX: custom( command )
% 
% Example: custom( 'print "hello world"' )
%
% See Also: http://lammps.sandia.gov/doc/Manual.html,
% http://lammps.sandia.gov/doc/Section_commands.html#comm

runCmd = InputFileElement();
runCmd.createInputFileText = @(~) customCommandCfg(command);

end

