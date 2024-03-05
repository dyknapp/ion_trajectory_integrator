function writeOutputStreamConfig( fHandle )
%WRITEOUTPUTSTREAMCONFIG

% Configures the output stream of the LAMMPS simulation to output every
% number of steps. This is a temporary workaround for a bug on windows
% machines where LAMMPS will not flush the buffer, causing diagnostic
% information to only be displayed at the end of the simulation for those
% machines. Intermittently writing text to the output buffer ensures that
% the output will be flushed.

fprintf(fHandle, '# Configuring additional output to flush buffer during simulation \n');
fprintf(fHandle, 'thermo 10000\n');
fprintf(fHandle, 'thermo_style custom step cpu\n');
fprintf(fHandle, 'thermo_modify flush yes\n');
fprintf(fHandle, 'fix extra all print 10 "SIMULATION RUNNING----------------------"\n\n');

end

