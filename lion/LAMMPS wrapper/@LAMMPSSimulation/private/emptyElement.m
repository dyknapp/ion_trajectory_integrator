function [ e ] = emptyElement(  )
%EMPTYELEMENT Creates an empty input file element

e = InputFileElement();
e.createInputFileText = @(~) {};
end

