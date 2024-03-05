function [ id ] = idString( n )
%IDSTRING Generates a random id string.

if nargin < 1
    n = 6;
end

    function c = randomChar()
        i = ceil(rand(1) * 26);
        i = min(i, 26); i = max(1,i);
        chars = 'a':'z';
        c = chars(i);
    end

id = [];
for j=1:n
    id = [id randomChar()];
end

end

