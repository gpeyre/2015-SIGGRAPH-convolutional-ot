function [newmesh indicator] = addCylinders(mesh,vf,vfBase,radius,nRadialSamples)

if nargin < 4
    radius = .0025;
end

if nargin < 5
    nRadialSamples = 5;
end

[X,T,indicator] = addCylindersHelper(mesh.vertices, mesh.triangles, vf, vfBase, radius, nRadialSamples);

newmesh.vertices = X;
newmesh.triangles = T;
newmesh.numVertices = size(X,1);