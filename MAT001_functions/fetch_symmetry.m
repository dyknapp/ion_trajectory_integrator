% A quick function for checking the symmetry of a SIMION file.
%
% dknapp, 01.11.2023

function symmetry = fetch_symmetry(path, start_line)
    fid = fopen(path,'r');
    for i = 1:start_line
     line = fgetl(fid);
     if contains(line, "symmetry")
         symmetry = extractAfter(line, "symmetry ");
         return
     end
    end
    fclose(fid);
end