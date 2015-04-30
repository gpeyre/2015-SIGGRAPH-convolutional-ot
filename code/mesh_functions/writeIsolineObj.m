function writeIsolineObj(M,data,filename,nLines,extra,radius)

if nargin < 5
    extra = .2;
end

if nargin < 5
    radius = .01; % need a better default...
end

rescaled = (data-min(data))/(max(data)-min(data));
[lineStart,lineEnd,index] = isolines(M.vertices,M.triangles,rescaled,linspace(0,1,nLines));

% showDescriptor(M,rescaled);
% hold on;
% plot3([lineStart(:,1) lineEnd(:,1)]',[lineStart(:,2) lineEnd(:,2)]',[lineStart(:,3) lineEnd(:,3)]','k-','linewidth',3);
% axis equal;

diff = (lineEnd-lineStart);

[M2 indicator] = addCylinders(M,(1+extra)*diff,lineStart-extra*diff/2,radius,8);
indicator = indicator*2;
indicator(1:M.numVertices) = rescaled;
showDescriptor(M2,indicator);

writeTexturedObj(filename, M2, indicator, 'level_set.mtl');