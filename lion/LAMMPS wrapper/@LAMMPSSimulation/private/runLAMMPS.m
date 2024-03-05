function [result] = runLAMMPS( inputFileName )
%RUNLAMMPS Runs lammps with the given file as an input argument.
result = false;

c = onCleanup(@() killProc());
    function killProc()
        if exist('process', 'var') && process.isAlive
            disp('Terminating lammps process.');
            process.destroy();
        end
    end

%Whatever error is thrown (or if user terminates), kill the lammps process
%if it has started.

[directoryS,~,~] = fileparts(mfilename('fullpath'));
[directoryS,~,~] = fileparts(directoryS);
[directoryS,~,~] = fileparts(directoryS);

%Try to find Lammps.cfg file in the wrapper folder.
if exist(fullfile(directoryS, 'Lammps.cfg'), 'file')
    fHandle = fopen(fullfile(directoryS, 'Lammps.cfg'), 'r');
    
    line = fgets(fHandle);
    while (line(1) == '#')
        line=fgets(fHandle);
    end
    lammpsPath = line;
    fclose(fHandle);
    
    command_string = "";
    command_string = command_string + lammpsPath + " ";
    command_string = command_string + "-in ";
    command_string = command_string + which(inputFileName);
    % command_string = command_string + "-np ";
    % command_string = command_string + sprintf('%s ', getenv('NUMBER_OF_PROCESSORS'));
    fprintf("%s\n", command_string);
    result = system(command_string)
    
    if result ~= 0
        error('LAMMPS executable did not run. Please ensure the program will run without errors in a terminal.');
    end
    result = true;
    
else
    error(['Cannot find Lammps.cfg file. Please make sure this exists'...
        ' within the folder `LAMMPS Wrapper`. This file should be a'...
        ' plain-text file containing the path to the lammps executable'...
        ' that LIon should run (eg C:\lammps\lmp_win_no-mpi.exe)']);
end

end

