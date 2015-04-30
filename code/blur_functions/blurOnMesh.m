function result = blurOnMesh(signal,M,time,steps,transpose)

h = time/steps;

% area*(f1-f0)/h = laplace*f1
% --> (area - h*laplace)f1 = area*f0

nv = M.numVertices;

% We really should pre-factor this for speed...
blurInverse = spdiags(M.areaWeights,0,nv,nv) - h*M.cotLaplacian;

result = signal;

if nargin < 5 || transpose == 0
    for i=1:steps
        result = blurInverse \ bsxfun(@times,result,M.areaWeights);
    end
else
    for i=1:steps
        result = bsxfun(@times,blurInverse' \ result, M.areaWeights);
    end
end