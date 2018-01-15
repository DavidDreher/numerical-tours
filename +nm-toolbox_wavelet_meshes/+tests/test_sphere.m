% test for spherical data processing

j = 5;
options.keep_subdivision = 0;
options.base_mesh = 'ico';
options.relaxation = 0;
[vertex,face] = nm-toolbox_wavelet_meshes.compute_semiregular_sphere(j,options);

options.relaxation = 4;
[vertex1,face] = nm-toolbox_wavelet_meshes.compute_semiregular_sphere(j,options);


clf;
subplot(1,2,1);
nm-toolbox_graph.plot_mesh(vertex,face);
shading faceted;
subplot(1,2,2);
nm-toolbox_graph.plot_mesh(vertex1,face);
shading faceted;
