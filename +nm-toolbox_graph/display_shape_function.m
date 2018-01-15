function B = display_shape_function(A,options)

% display_shape_function - display a function inside a shape.
%
%   B = display_shape_function(A,options);
%
%   options.cm gives colormap.
%   options.display_levelsets = 1 to show level sets.
%   options.pstart to give location of start point.
%
%   Copyright (c) Gabriel Peyre 2010

options.null = 0;

cm = nm-toolbox_general.getoptions(options, 'cm', jet(256));
display_levelsets = nm-toolbox_general.getoptions(options, 'display_levelsets', 0);
nbr_levelsets = nm-toolbox_general.getoptions(options, 'nbr_levelsets', 20);
pstart = nm-toolbox_general.getoptions(options, 'pstart', []);

A(isinf(A)) = 0;
I = find(A>0);
A(I) = nm-toolbox_general.rescale(A(I));

n = size(A,1);
B = ones(n,n,3);

m = size(cm,1);


for i=1:3
    c = cm(:,i);
    Bi = B(:,:,i);   
    Bi(I) = c( round(A(I)*(m-1))+1 );
    B(:,:,i) = Bi;
end

if nargout==0
	hold on;
    nm-toolbox_signal.imageplot(B);
    if display_levelsets
        contour(nm-toolbox_general.rescale(A), nbr_levelsets, 'k');
        % set(gca, 'LineWidth', 2);
    end
    if not(isempty(pstart))
        h = plot(pstart(2), pstart(1), 'r.');
        set(h, 'MarkerSize', 30);
    end
    hold off;
    axis ij;
end
