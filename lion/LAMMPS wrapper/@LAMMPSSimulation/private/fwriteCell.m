function fwriteCell( f, cell )
%FWRITECELL writes a cell of strings to a file as lines.

for i=1:length(cell)
   fprintf(f,'%s\n', cell{i}); 
end

end

