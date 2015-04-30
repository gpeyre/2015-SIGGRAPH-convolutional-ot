function [Xout,Tout,index] = addSpheres(X,T,whichVerts,radii)

if nargin < 4
    radii = ones(length(whichVerts),1);
end

[Xsphere,Tsphere] = readOff('simpleSphere.off');

Xout = X;
Tout = T;
index = zeros(size(X,1),1);

for i=1:length(whichVerts)
    location = X(whichVerts(i),:);
    shiftedX = Xsphere * radii(i) + repmat(location,size(Xsphere,1),1);
    
    Tout = [Tout;Tsphere + size(Xout,1)];
    Xout = [Xout;shiftedX];
    index = [index;i*ones(size(Xsphere,1),1)];
end