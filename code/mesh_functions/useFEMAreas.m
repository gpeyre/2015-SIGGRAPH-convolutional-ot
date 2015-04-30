function M = useFEMAreas(M,numEigs,FEM)

% eigenvalues should be negative

if numEigs ~= M.numVertices
    [evecs, evals] = eigs(FEM.laplacian,FEM.vtxInnerProds,max(numEigs,1), 1);
else
    fprintf('Using full basis.\n');
    [evecs, evals] = eig(full(FEM.laplacian),full(FEM.vtxInnerProds));
end

evals = diag(evals);
% norm(FEM.laplacian*evecs - bsxfun(@times,FEM.vtxInnerProds*evecs,evals'),'fro')

M.laplaceBasis = evecs;
M.eigenvalues = evals;
M.faceAreas = FEM.faceAreas;

% [minEig,idx] = min(abs(M.eigenvalues)); % make this the first column
% 
% tempEig = M.eigenvalues(1);
% M.eigenvalues(1) = minEig;
% M.eigenvalues(idx) = tempEig;
% 
% tempVec = M.laplaceBasis(:,idx);
% M.laplaceBasis(:,idx) = M.laplaceBasis(:,1);
% M.laplaceBasis(:,1) = tempVec;

[~,order] = sort(abs(M.eigenvalues));
M.eigenvalues = M.eigenvalues(order);
M.laplaceBasis = M.laplaceBasis(:,order);

norms = sum((FEM.vtxInnerProds*M.laplaceBasis).*M.laplaceBasis,1);
M.laplaceBasis = bsxfun(@rdivide,M.laplaceBasis,sqrt(norms));

% norm(FEM.laplacian*M.laplaceBasis-FEM.vtxInnerProds*bsxfun(@times,M.laplaceBasis,M.eigenvalues'),'fro')