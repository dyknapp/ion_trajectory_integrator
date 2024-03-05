function [ output ] = cfghelperTimestep( input )
%cfghelperTimestep Gets or sets the global variable used when writing the
%config file. An end user of the wrapper should never need to use this
%function.

%This function is a little odd and so warrants further discussion:
% The intended use case of this function is to allow implementation of
% variable timestep size into LIon so that a single simulation may have
% different timesteps as fixes are added/removed eg before/after
% minimisation.
% 
% To implement this, the cfghelperTimestep function is used, which allows
% access to a persistent variable. The function is only invoked when
% writing the input file. Whenever the timestep changes in lammps, the
% persistent variable is also changed. Functions requiring access to the
% current timestep size (eg for time averaging purposes, or to ensure we
% simulate for x milliseconds) can themselves invoke the 'get' of
% cfghelperTimestep when writing their input file text to ensure the
% current value of the timestep is used.
%
% As a result, this function probably isn't what you are looking for.
% Perhaps you mean to use 'sim.TimeStep' (current, until variable timestep
% introduced)?

persistent timestep;

if ischar(input) && strcmpi(input, 'get')
   output = timestep;
else
    timestep = input;
end

end

