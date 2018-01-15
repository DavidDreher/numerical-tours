function f = nm-toolbox_graph.perform_mesh_smoothing(face,vertex,f,options)

% perform_mesh_smoothing - smooth a function defined on a mesh by averaging
%
%   f = perform_mesh_smoothing(face,vertex,f,options);
%
%   Smooth a function f on a width of options.niter_averaging vertices.
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
naver = nm-toolbox_general.getoptions(options, 'niter_averaging', 1);
type = nm-toolbox_general.getoptions(options, 'averaging_type', 'combinatorial');

if nargin<3
    f = [];
end
if isempty(f)
    f = vertex;
end
if size(f,1)<size(f,2)
    f = f';
end
[vertex,face] = nm-toolbox_graph.check_face_vertex(vertex,face);

if size(f,2)>1
    for i=1:size(f,2)
        f(:,i) = nm-toolbox_wavelet_meshes.toolbox.nm-toolbox_graph.perform_mesh_smoothing(face,vertex,f(:,i),options);
    end
    return;
end

n = max(face(:));

% compute normalized averaging matrix
if strcmp(type, 'combinatorial')
    %add diagonal
    W = nm-toolbox_graph.triangulation2adjacency(face) + speye(n);
    D = spdiags(full(sum(W,2).^(-1)),0,n,n);
    W = D*W;
else
    options.normalize=1;
    W = nm-toolbox_graph.compute_mesh_weight(vertex,face,type,options);
end

% do averaging to smooth the field
for k=1:naver
    f = W*f;
end
