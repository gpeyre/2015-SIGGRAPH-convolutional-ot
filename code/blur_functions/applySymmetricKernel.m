function result = applySymmetricKernel(signal,M,time,steps)

h = time/steps;

nv = M.numVertices;

% We really should pre-factor this for speed...
blurInverse = spdiags(M.areaWeights,0,nv,nv) - h*M.cotLaplacian;
result = bsxfun(@rdivide,signal,M.areaWeights); % post-multiply kernel by 1/areas

for i=1:steps
    result = blurInverse \ bsxfun(@times,result,M.areaWeights);
end