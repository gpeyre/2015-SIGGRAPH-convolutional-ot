function generateSoftMapMeshes(Msource,Mtarget,sourceFile,targetFile,fixedSources,fixedTargets,sourceFunctionPoints,f,sphereRadius)

% First make source mesh

[X1,T1,idx2] = addSpheres(Msource.vertices,Msource.triangles,fixedSources,sphereRadius*ones(size(fixedSources)));
[X2,T2,idx] = addSpheres(X1,T1,sourceFunctionPoints,sphereRadius*ones(size(sourceFunctionPoints)));

M.vertices = X2;
M.triangles = T2;
M.numVertices = size(X2,1);

idx(idx2 ~= 0) = 8;
%idx(1:length(idx2)) = idx2;
showDescriptor(M,idx);

writeTexturedObj(sourceFile, M, idx, 'source.mtl');

% Now make target mesh

[X,T,idx] = addSpheres(Mtarget.vertices,Mtarget.triangles,fixedTargets,sphereRadius*ones(size(fixedTargets)));

f = bsxfun(@times,f,min(f));
f = bsxfun(@rdivide,f,max(f));
f = [zeros(size(f,1),1) f];
f((end+1):size(X,1),:) = 0;
f(idx ~= 0,9) = 2;

M.vertices = X;
M.triangles = T;
M.numVertices = size(X,1);

[~,test] = max(f,[],2);
showDescriptor(M,test)

writeTexturedObj(targetFile, M, f, 'softmap.mtl', 9);