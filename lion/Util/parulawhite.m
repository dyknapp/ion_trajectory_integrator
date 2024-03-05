function cmap = parulawhite
%PARULAWHITE Parula white colormap, with white for zeros.

cmap = parula(8);
cmap(1,:) = [ 0.9 0.9 0.9 ];
cmap(2,:) = [ 0.8 0.8 0.8 ];
cmap = interp1(linspace(0, 1, 8), cmap, linspace(0,1,256));

end