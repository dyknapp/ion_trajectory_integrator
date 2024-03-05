function [ strings ] = outputCfg( ID, filename, outputFreq, style, args )
%OUTPUTCFG Cfg file content for outputting data from lammps.

strings = { sprintf('dump %d all %s %d %s %s', ID, style, outputFreq, filename, args)};

end

