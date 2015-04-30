function FEM = firstOrderFEM(mesh,basisSize)

if nargin < 2
    basisSize = 0;
end

% This gives us a struct with fields
% 
% FEM.faceAreas
% FEM.vtxInnerProds
% FEM.faceInnerProds
% FEM.grad{i} for i = 1,2,3
% FEM.laplacian
% 
% and optionally, the fields
%
% FEM.laplaceEigenvalue
% FEM.laplaceBasis

% Read off per-triangle vertex locations
for i=1:3
    faceVtx{i} = mesh.vertices(mesh.triangles(:,i),:);
end

% Compute edge lengths for the triangles
% (some calculations repeated here, whatever)
for i=1:3
    j = mod(i,3)+1;
    d = faceVtx{j} - faceVtx{i};
    edgeLengths{i} = sqrt(sum(d.^2,2));
end

% Compute face areas using Heron's formula
% (again, not the fastest/most stable strategy)
s = (edgeLengths{1}+edgeLengths{2}+edgeLengths{3})/2;
faceAreas = sqrt(s.*(s-edgeLengths{1}).*(s-edgeLengths{2}).*(s-edgeLengths{3}));
%totalArea = sum(faceAreas);
%faceAreas = faceAreas / totalArea;
FEM.faceAreas = faceAreas;

% edges(:,1) < edges(:,2)
r = [];
c = [];
v = [];
for i = 1:3
    j = mod(i,3)+1;
    
    % Inner product of any two hat functions is area/12 on a given face
    t1 = mesh.triangles(:,i);
    t2 = mesh.triangles(:,j);
    r = [r; t1; t2];
    c = [c; t2; t1];
    v = [v; faceAreas/12; faceAreas/12];
    
    % Self product is area/6
    r = [r; t1];
    c = [c; t1];
    v = [v; faceAreas/6];
end

FEM.vtxInnerProds = sparse(r,c,v);
FEM.faceInnerProds = sparse(1:mesh.numTriangles, 1:mesh.numTriangles, faceAreas);

% We'll make three matrices representing the x,y,z components of the FEM.grad
N = cross(faceVtx{3}-faceVtx{1},faceVtx{2}-faceVtx{1});
N = bsxfun(@rdivide,N,sqrt(sum(N.^2,2)));
r = {[],[],[]};
c = {[],[],[]};
v = {[],[],[]};
for i=1:3
    j = mod(i,3)+1;
    k = mod(j,3)+1;
    
    d = faceVtx{k} - faceVtx{j}; % check the sign!
    rot = cross(d,N);
    scaled = bsxfun(@rdivide,rot,2*faceAreas); %is this scaling right? -- used to divide rot by /sqrt(totalArea)
    
    for j = 1:3
        r{j} = [r{j} ; (1:mesh.numTriangles)'];
        c{j} = [c{j} ; mesh.triangles(:,i)];
        v{j} = [v{j} ; scaled(:,j)];
    end
end

FEM.grad = {};
for i = 1:3
    FEM.grad{i} = sparse(r{i},c{i},v{i});
end

FEM.laplacian = -FEM.grad{1}'*FEM.faceInnerProds*FEM.grad{1} - FEM.grad{2}'*FEM.faceInnerProds*FEM.grad{2} - FEM.grad{3}'*FEM.faceInnerProds*FEM.grad{3};

% Make rotated gradient -- N cross grad
FEM.rotGrad{1} = bsxfun(@times,mesh.faceNormals(:,2),FEM.grad{3}) - bsxfun(@times,mesh.faceNormals(:,3),FEM.grad{2});
FEM.rotGrad{2} = bsxfun(@times,mesh.faceNormals(:,3),FEM.grad{1}) - bsxfun(@times,mesh.faceNormals(:,1),FEM.grad{3});
FEM.rotGrad{3} = bsxfun(@times,mesh.faceNormals(:,1),FEM.grad{2}) - bsxfun(@times,mesh.faceNormals(:,2),FEM.grad{1});

if basisSize >= 1
    [B,D] = eigs(FEM.laplacian, basisSize, 'sm');
    integrals = diag(B'*FEM.vtxInnerProds*B);
    FEM.laplaceEigenvalues = diag(D);
    FEM.laplaceBasis = bsxfun(@rdivide,B,sqrt(integrals)');
end