%%
% Test for OT color transfer.

addpath('../toolbox/');
addpath('../colors_functions/');
addpath('../image_blur/');
addpath('../blur_functions//');
addpath('../convolutional_wasserstein/');
addpath('../../data/meshes/');
addpath('../mesh_functions/');

% #bins
N = 40;
N = 60;

%%
% helpers

mmin = @(x)min(x(:));
mmax = @(x)max(x(:));
normalize = @(h)h/sum(h(:));
setfigname = @(name)set(gcf, 'Name', name, 'NumberTitle','off');
% plot histograms
delta = .1/N^2;
dispHist = @(x)-log(x+delta);
dispCell = @(H)cellfun(dispHist, H, 'UniformOutput', false);
Entropy = @(x)-sum(x(x>0).*log(x(x>0)));

%%
% Blur kernel

mu = N/50;
mu = N/40;
blur = load_filtering('imgaussian', N);
K = @(x)blur(x,mu);
Kv = @(x)apply_3d_func(K,x);

%%
% Load densities

if 1 % not(exist('names'))
    %names = {'spiky', 'spheres'};
    %names = {'spiky', 'sphere'};
    %names = {'sphere', 'boxes'};
    %names = {'duck', 'horse'};
    %names = {'duck' 'horse' 'shears' 'moomoo_s0'};
    names = {'spiky','sphere','boxes','cone_rotated'};
end
p = length(names);

f = {};
for i=1:p
    name = names{i};
    f{i} = normalize( load_volume(names{i}, N) );
    if isempty(f{i})
        %% try to load a mesh %%
        [V,F] = read_off([name '.off']);
        V = rescale(V, .03, .97);
        t = linspace(0,1,N);
        f{i} = inpolyhedron(F',V',t,t,t);
    end
end


opts.alpha = 1; % transparency
opts.color = [0 1 0];
clf;
for i=1:p
    subplot(1,length(names),i);
    plot_isosurface(f{i},opts);
end

rep = ['../results/volumetric/' names{1}];
for i=2:p
    rep = [rep '-' names{i}];
end
rep = [rep '/'];

if not(exist(rep))
    mkdir(rep);
end


%% 
% Compute transport coupling
%   pi = diag(w1)*K*diag(w0)  
% and  
%   pi*1 = w1.*K(w0) = p1 
% and 
%   pi'*1 = w1.*K(w0) = p0

options.tol = 1e-9;
options.tol = 0;
options.niter = 100;
options.verb = 2;
slicing = @(x)x(:,:,end/2);
options.disp = @(w0,w1)imageplot( slicing(w0.*K(w1)) );

% clf; 
% [distances,w0,w1] = convolutionalDistance(f{1}, f{2}, [], K,[], options);

%%
% Compute displacemet interpolation

q = 5; 
switch p
    case 2
        % displacement interpolation
        t = linspace(0,1,q);
        W = [t;1-t];
    case 3
        % triangle interpolation
        W = [ ...
            [0, 0, 1]; ...
            [1, 0, 3]; [0, 1, 3]; ...
            [1,0,1]; [1,1,2]; [0,1,1]; ...
            [3,0,1]; [2,1,1]; [1,2,1]; [0,3,1]; ...
            [1,0,0]; [3,1,0]; [1,1,0]; [1,3,0]; [0,1,0] ...
            ]';
    case 4
        % bilinear interpolation
        t = linspace(0,1,q);
        [T,S] = meshgrid(t,t); S = S(:); T = T(:);
        W = [(1-S).*(1-T) S.*(1-T) (1-S).*T S.*T]';
end
Q = size(W,2);

cachedMeshes = cell(Q,1);

parfor i=1:Q
    fprintf('%d of %d...\n',i,Q);
    w = W(:,i)'; w = w/sum(w);
    
    % select entropy bound
	entropyLimit = [];
    % store as 2D matrix
    Hv = []; 
    for k = 1:p        
        Hv = [Hv f{k}(:)];
    end
    % do the computation
    options = [];
    options.disp = @(x) plot_isosurface(reshape(x,[N N N]));
    options.verb = 2;
    options.tol = 0;
    options.niter = 100;
    [B,u] = convolutionalBarycenter( Hv, w, [], Kv, [],entropyLimit, options);
    B = reshape(B, [N N N]);
    B = B/max(B(:));
    % display
%     clf; 
    opts = [];
    opts.alpha = 1; % transparency
    if p==2
        opts.color = [w(1) 0 w(2)];
    else
        opts.color = [w(1) w(2) w(3)];
    end
    opts.isolevel = median(B(:));
    opts.isolevel = (max(B(:))-min(B(:)))/2;
    F = plot_isosurface(B,opts);    
    close all;
    % save as image
%     if p==4
%         str = ['barycenter-' num2str(S(i)*(q-1)) '-' num2str(T(i)*(q-1))];
%     else
%         str = ['barycenter-' num2str(i)];
%     end
%     saveas(gcf, [rep str '.png'], 'png');

    % save as object
    mesh = [];
    mesh.vertices = F.vertices;
    mesh.triangles = F.faces;
    mesh.numVertices = size(F.vertices, 1);
%     vertexTexture = zeros(mesh.numVertices,1);
%     writeTexturedObj([rep 'barycenter-' str '.obj'], mesh, vertexTexture); 
    
    cachedMeshes{i} = mesh;
    
%     shift = [0 170*S(mod(i,q)+1) 170*T(floor(i/q+1))];
%     fullX = [fullX ; bsxfun(@plus,mesh.vertices,shift)];
%     if i>1
%         fullT = [fullT ; mesh.triangles+max(fullT(:))];
%     else
%         fullT = mesh.triangles;
%     end
end
 
%% Rescale to unit box

minCoord = Inf;
for i=1:length(cachedMeshes)
    minCoord = min(minCoord,min(cachedMeshes{i}.vertices(:)));
end

maxCoord = -Inf;
for i=1:length(cachedMeshes)
    cachedMeshes{i}.vertices = cachedMeshes{i}.vertices-minCoord;
    maxCoord = max(maxCoord,max(cachedMeshes{i}.vertices(:)));
end

for i=1:length(cachedMeshes)
    cachedMeshes{i}.vertices = cachedMeshes{i}.vertices/maxCoord;
end

%% Combine meshes

fullX = [];
fullT = [];
textureCoords = [];

for i=1:length(cachedMeshes);
    fullX = [fullX ; bsxfun(@plus,cachedMeshes{i}.vertices,.85*[S(i)*q,T(i)*q,0])];
    if i>1
        fullT = [fullT ; cachedMeshes{i}.triangles+max(fullT(:))];
    else
        fullT = cachedMeshes{i}.triangles;
    end
    
    textureCoords = [textureCoords; repmat([S(i) T(i)],cachedMeshes{i}.numVertices,1)];
end

MM.vertices = fullX;
MM.triangles = fullT;
MM.numVertices = size(fullX,1);

%% Make .obj with constant color

filename = 'gridinterp.obj';
materialFile = '2dcolor.mtl';
fid = fopen(filename,'wt');

fprintf(fid, 'mtllib %s\n', materialFile);
fprintf(fid, 'usemtl material_0\n');

fprintf(fid, 'v %f %f %f\n', MM.vertices');
fprintf(fid, 'vt %f %f\n', .05+.9*textureCoords');

face = MM.triangles;
face_texcorrd = [face(:,1), face(:,1), face(:,2), face(:,2), face(:,3), face(:,3)];
fprintf(fid, 'f %d/%d %d/%d %d/%d\n', face_texcorrd');

fclose(fid);