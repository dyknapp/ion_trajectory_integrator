function plot_img(xs, ys, image)
    imagesc(xs, ys, image);
    axis image; set(gca,'YDir','normal'); colormap('pink');
end