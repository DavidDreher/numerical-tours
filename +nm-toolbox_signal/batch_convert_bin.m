% convert image to .bin
name = 'flowers';
M = nm-toolbox_signal.load_image(name);
nm-toolbox_signal.write_bin(nm-toolbox_general.rescale(M), name);
