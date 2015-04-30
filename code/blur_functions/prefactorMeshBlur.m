function structure = prefactorMeshBlur(M,time,steps,batchsize)

if nargin < 4
    batchsize = 500;
end

structure.batchsize = batchsize;

h = time/steps;

% area*(f1-f0)/h = laplace*f1
% --> (area - h*laplace)f1 = area*f0
% --> (area - h*laplace)(f1-f0) = h*laplace*f0
% --> (area/h - laplace)(f1-f0) = laplace*f0 

% I use the last one -- the matrix is worse conditioned but I think should
% capture f1-f0 in more detail, which is important for small h.

nv = M.numVertices;

% We really should pre-factor this for speed...
structure.blurInverse = spdiags(M.areaWeights,0,nv,nv)/h - M.cotLaplacian;%spdiags(M.areaWeights,0,nv,nv) - h*M.cotLaplacian;
structure.ordering = symamd(structure.blurInverse);
structure.R = chol(structure.blurInverse(structure.ordering,structure.ordering));
structure.steps = steps;
structure.reorderedAreaWeights = M.areaWeights(structure.ordering);
structure.reorderedLaplacian = M.cotLaplacian(structure.ordering,structure.ordering);