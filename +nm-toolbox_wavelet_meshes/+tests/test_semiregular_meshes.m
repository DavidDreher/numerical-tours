% test for semi-regular meshes

path(path, 'toolbox/');
path(path, 'gim/');

name = 'bunny';
name = 'cow';
name = 'armadillo';
name = 'horse';
name = 'gargoyle';

options.func = 'mesh';
options.name = name;
options.use_elevation = 0;
options.use_color = 0;

%% Load the geometry image
M = nm-toolbox_wavelet_meshes.read_gim([name '-sph.gim']);
%% create the semi regular mesh from the Sph-GIM
J = 7;
[vertex,face,vertex0] = nm-toolbox_wavelet_meshes.compute_semiregular_gim(M,J,options);


rep  = 'results/semi-regular/';
if not(exist(rep))
    mkdir(rep);
end

for j=3:J
    clf;
    nm-toolbox_wavelet_meshes.plot_spherical_function(vertex{j},face{j},[], options);
    shading faceted;
    pause(1);
    saveas(gcf, [rep name '-scale-' num2str(j) '.png'], 'png');
end
