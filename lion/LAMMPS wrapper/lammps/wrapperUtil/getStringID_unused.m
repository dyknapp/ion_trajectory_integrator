function [ string ] = getStringID( id )
%GETSTRINGID Converts an integer ID into a string of characters.
% This converts the given number into a representation of base 'a':'z'. It
% is used to write IDs in Lammps when identifiers must be string
% characters, not numbers.
% 
%Example:
%  getStringID(12345)

symbols = 'a':'z';
string = dec2base(id, length(symbols));

%convert numbers to string ids
for i=1:length(string)
   c = string(i);
   if (c < 'A')
       c = c + 'A'-'0';
   else
       c = c+10;
   end
   string(i) = char(c);
end

end