function [Blur, Delta,Grad] = blurAnisotropic(M, options)

% blurAnisotropic - load anisotropic heat kernel
%
%   [Blur, Delta, Grad] = blurAnisotropic(H, t, filtIter, options);
%
%   H is an (N,N,2,2) representing the 2D Riemannian manifold.
%
%   Blur(u,t,filtIter) iterated filtIter times
%       u = (Id + t*(Delta))\u
%
%   Copyright (c) 2014 Gabriel Peyre

if nargin<2
    filtIter = 1;
end

N = size(M,1);

id = speye(N);

options.null = 0;
diff_type = getoptions(options, 'diff_type', 'fwd');
e = ones(N,1);
switch diff_type
    case 'sym'
        % symmetric
        D = spdiags([e -e]/2, [-1 1], N, N);
        D([1 end],:) = 0; % Neumann BC
    case 'fwd'
        % forward
        D = spdiags([e -e], [0 1], N, N);
        D(end,:) = 0; % Neumann BC
    otherwise
        error('Unknown finite differenciation mode.');
end

Dx = kron(D,id);
Dy = kron(id,D);
Grad = [Dx;Dy];

% metric along the diagonal
diagH = @(h)spdiags(h(:), 0, N^2,N^2);
dH = [diagH(M(:,:,2,2)), diagH(M(:,:,1,2)); ...
      diagH(M(:,:,2,1)), diagH(M(:,:,1,1))];
% Laplacian
Delta = Grad'*dH*Grad;
% load blurring kernel
Blur = @(u,t,filtIter)heat_iter(Delta, t, filtIter, u);

end

%%%
function u = heat_iter(Delta, t, filtIter, u)

s = size(u);
u = u(:);
Id = speye(size(u,1));
for i=1:filtIter
    u = (Id + t*Delta)\u;
end
u = reshape(u,s);

end
