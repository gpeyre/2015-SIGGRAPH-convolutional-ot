function M = fit_metric(fDir)

% fit_metric - fit a 2D Riemannian metric
%
%   M = fit_metric(fDir);
%
%   fDir(h) gives the value of the metric in direction h (a 1x2 vector). 
%   M is a 2x2 symmetric matrix, that fit the measurements, i.e. such that
%       fDir(h)^2 ~ h*M*h'
%
%   Copyright (c) 2014 Gabriel Peyre



a = fDir([1 0])^2;
b = fDir([0 1])^2;
c = (fDir([1 1])^2-a-b)/2;
M = [a c; c b];

[U,S] = eig(M); S = diag(S);
if min(S)<0
    warning('Non positive matrix, clamping singular values');
    M = U*diag(max(S,0))*U';
end