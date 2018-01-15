function y = callback_optical_flow(x,tmp)

% callback_optical_flow - callback for the optical flow regularization.
%
%   y = callback_optical_flow(x);
%
%   This is a symmetric definite positive operator.
%
%   Copyright (c) 2009 Gabriel Peyre

global lambda; global D;
n = size(D,1);
v = reshape(x, [n n 2]);
u = sum(v.*D,3);
y = D.*repmat(u, [1 1 2]) - lambda * cat(3, nm-toolbox_signal.div(nm-toolbox_signal.grad(v(:,:,1))), nm-toolbox_signal.div(nm-toolbox_signal.grad(v(:,:,2))));
y = y(:);
