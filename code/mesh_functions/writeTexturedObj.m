function writeTexturedObj(filename, mesh, vertexTexture, materialFile, totalTextures)

% overide default displaying
do_display = 0; 

if nargin < 4
    materialFile = 'scalarFunction.mtl';
end

if nargin < 5
    totalTextures = 1;
end

vertex = mesh.vertices;
face = mesh.triangles;

if (size(vertex,2)~=3), vertex=vertex'; end
if size(vertex,2)~=3, error('vertex does not have the correct format.'); end

if size(face,2)~=3, face=face'; end
if size(face,2)~=3, error('face does not have the correct format.'); end

fid = fopen(filename,'wt');
if( fid==-1 ), error('Can''t open the file.'); end

fprintf(fid, 'mtllib %s\n', materialFile);
fprintf(fid, 'usemtl material_0\n');

% vertex position
fprintf(fid, 'v %f %f %f\n', vertex');

top = max(vertexTexture,[],1);
bottom = min(vertexTexture,[],1);
diff = top-bottom;
vertexTexture = (vertexTexture - repmat(bottom,mesh.numVertices,1))./repmat(diff,mesh.numVertices,1); % map to [0,1]
vertexTexture = vertexTexture * .99 + .005;

[maxTexture,index] = max(vertexTexture,[],2);

spacing = 1/totalTextures;

% ones(mesh.numVertices,1);
% for i=1:totalTextures
%     mask = double(vertexTexture(:,i) == maxTexture);
%     index = index.*(1-mask) + i * mask;
% end

ycoord = (index*spacing - spacing/2);
fprintf(fid, 'vt %f %f\n', [maxTexture ycoord]');

if do_display
    showDescriptor(mesh,ycoord);
    showDescriptor(mesh,maxTexture);
end

% face
face_texcorrd = [face(:,1), face(:,1), face(:,2), face(:,2), face(:,3), face(:,3)];
fprintf(fid, 'f %d/%d %d/%d %d/%d\n', face_texcorrd');

fclose(fid);
