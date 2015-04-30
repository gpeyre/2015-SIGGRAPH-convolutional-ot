function h = showDescriptor(mesh, f, title, camera, scale, h)
if nargin < 3
    title = '';
end

if nargin < 4
    camera = [];
end

if nargin < 5
    scale = [];
end

if nargin < 6
    h = figure;
end

X = mesh.vertices;
T = mesh.triangles;
parentFig = []; %getParentFigure(h);

isAxes = strcmp('axes',get(h,'type'));

if isAxes
    if size(f,1) == mesh.numVertices
        patch('vertices',X,'Faces',T,'FaceColor','interp','CData',double(f),'edgecolor','none','parent',h); 
    else
        patch('vertices',X,'Faces',T,'FaceColor','flat','CData',double(f),'edgecolor','none','parent',h); 
    end
    axis(h,'equal');
    colorbar('peer',h);
    cameratoolbar;%(parentFig,'Show');
else % h is a figure
    figure(h);
    if size(f,1) == mesh.numVertices
        patch('vertices',X,'Faces',T,'FaceColor','interp','CData',double(f),'edgecolor','none'); 
    else
        patch('vertices',X,'Faces',T,'FaceColor','flat','CData',double(f),'edgecolor','none'); 
    end
    axis equal;
    colorbar;
    cameratoolbar;
end

if ~isempty(title)
    set(h,'name',title);
end

if nargin >= 4 && isfield(camera, 'up')
    campos(camera.position);
    camva(camera.angle);
    camup(camera.up);
end

if nargin >= 5 && length(scale) == 2
    caxis(scale);
end

if isAxes
    axis(h,'off');
else
    axis off;
end
