function coords = meanValueCoordinates(mesh, x)
% Implementation of Mean Value Coordinates for Closed Triangular Meshes
% Tao Ju, Scott Schaefer, Joe Warren

tol = 1e-4;

nv = mesh.numVertices;
nf = mesh.numTriangles;

d = normv(bsxfun(@minus, mesh.vertices, x));
d1 = d(mesh.triangles(:,1),:);
d2 = d(mesh.triangles(:,2),:);
d3 = d(mesh.triangles(:,3),:);

[minVal,minPos] = min(d);

if minVal < tol
    coords = zeros(nv,1);
    coords(minPos) = 1;
    return;
end

u = bsxfun(@rdivide,bsxfun(@minus, mesh.vertices, x),d);
u1 = u(mesh.triangles(:,1),:);
u2 = u(mesh.triangles(:,2),:);
u3 = u(mesh.triangles(:,3),:);

l = zeros(nf,3);
l(:,1) = normv(u2-u3);
l(:,2) = normv(u1-u3);
l(:,3) = normv(u1-u2);

theta = 2*asin(l/2);
h = sum(theta,2)/2;

specialW = zeros(nf,3);
specialW(:,1) = sin(theta(:,1)).*d2.*d3;
specialW(:,2) = sin(theta(:,2)).*d1.*d3;
specialW(:,3) = sin(theta(:,3)).*d2.*d1;

c(:,1) = 2*sin(h).*sin(h-theta(:,1)) ./ (sin(theta(:,2)).*sin(theta(:,3))) - 1;
c(:,2) = 2*sin(h).*sin(h-theta(:,2)) ./ (sin(theta(:,1)).*sin(theta(:,3))) - 1;
c(:,3) = 2*sin(h).*sin(h-theta(:,3)) ./ (sin(theta(:,2)).*sin(theta(:,1))) - 1;

det  = u1(:,1).*u2(:,2).*u3(:,3) + u1(:,2).*u2(:,3).*u3(:,1) + u1(:,3).*u2(:,1).*u3(:,2) - u1(:,3).*u2(:,2).*u3(:,1) - u1(:,1).*u2(:,3).*u3(:,2) - u1(:,2).*u2(:,1).*u3(:,3);
sign_det = sign(det);
s = bsxfun(@times,sign_det,sqrt(1 - c.^2));

w = zeros(nf,3);
w(:,1) = (theta(:,1) - c(:,2).*theta(:,3) - c(:,3).*theta(:,2)) ./ (d1.*sin(theta(:,2)).*s(:,3));
w(:,2) = (theta(:,2) - c(:,3).*theta(:,1) - c(:,1).*theta(:,3)) ./ (d2.*sin(theta(:,3)).*s(:,1));
w(:,3) = (theta(:,3) - c(:,1).*theta(:,2) - c(:,2).*theta(:,1)) ./ (d3.*sin(theta(:,1)).*s(:,2));

unacceptableS = sum(abs(s)<tol,2)~=0;
w(unacceptableS,:) = 0;

useSpecialW = find(pi-h<tol);
if ~isempty(useSpecialW)
    %w(useSpecialW,:) = specialW(useSpecialW,:);
    w = zeros(nf,3);
    w(useSpecialW(1),:) = specialW(useSpecialW(1),:);
end

weights = accumarray(mesh.triangles(:),w(:));
 
coords = weights ./ sum(weights);

function n = normv(v)
n = sqrt(sum(v.^2,2));