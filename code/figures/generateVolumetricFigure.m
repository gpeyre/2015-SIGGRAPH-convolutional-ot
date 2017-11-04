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
mu = N/25;
blur = load_filtering('imgaussian', N);
K = @(x)blur(x,mu);
Kv = @(x)apply_3d_func(K,x);

%%
% Load densities

if 1 % not(exist('names'))
    names = {'spiky', 'spheres'};
    names = {'spiky', 'sphere'};
    names = {'sphere', 'boxes'};
    names = {'duck', 'horse'};
    names = {'duck' 'horse' 'shears' 'moomoo_s0'};
    names = {'spot_flipped','duck','torus'};
    names = {'spiky','sphere','boxes','cone_rotated'};
    names = {'duck','spiky','moomoo_s0','double-torus'};
    names = {'shears' 'chair' 'wolf0' 'duck' };
    names = {'mushroom' 'torus' 'hand1' 'trim-star'}; % 'star_subdivided'
    names = {'torus'  'duck' 'mushroom' 'spiky'}; % 'star_subdivided'
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
        % check if rotation is needed
        A = eye(3);
        switch name
            case 'duck'
                A = [-0.0193 -0.2900 0.9568; 0.9677 0.2352 0.0908; -0.2514  0.9277 0.2761];
            case 'moomoo_s0'
                A = [-0.0971 0.7483 -0.6562; -0.9928 -0.1191 0.0111; -0.0699 0.6525 0.7545];
            case 'homer'
                A = [-0.4287  0.5566 -0.7116; -0.7577 -0.6505 -0.0524; -0.4921 0.5167 0.7006];
            case 'dinosaur'
                A = [-0.8362    0.0082   -0.5483   -0.2390   -0.9054    0.3509   -0.4935    0.4245    0.7591];
            case 'chair'
                A = [-0.8165    0.2728   -0.5089    0.3447   -0.4768   -0.8086   -0.4632   -0.8356    0.2952];
            case 'trim-star'
                A = [-0.7297    0.3243    0.6020   -0.6814   -0.4187   -0.6003    0.0574   -0.8483    0.5265];
            case 'hand1'
                A = [-0.4030   -0.7111   -0.5761    0.8945   -0.4391   -0.0837   -0.1934   -0.5491    0.8131];
            case 'mushroom'
                A = [-0.1212   -0.9430    0.3099    0.9227    0.0081    0.3855   -0.3660    0.3326    0.8691];
        end
        V = reshape(A,[3 3])*(V-.5)+.5;
        %  V0 = V; [A,~] = qr(randn(3)); V = A*(V0-.5)+.5;
        t = linspace(0,1,N);
        f{i} = inpolyhedron(F',V',t,t,t);
        opts.isolevel = .5;
        clf; plot_isosurface(f{i},opts);
        a=1;
    end
end


opts.alpha = 1; % transparency
opts.color = [0 1 0];
opts.isolevel = .5;
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
% Compute displacement interpolation

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

%%
% Do the computations.

for i=1:Q
    w = W(:,i)'; w = w/sum(w);
    % select entropy bound
	entropyLimit = [];
    % store as 2D matrix
    Hv = []; 
    for k = 1:p        
        Hv = [Hv f{k}(:)/sum(f{k}(:))];
    end
    % do the computation
    options.disp = @(x)plot_isosurface(reshape(x,[N N N]));
    [B,u] = convolutionalBarycenter( Hv, w, [], Kv, [],entropyLimit, options);
    B = reshape(B, [N N N]);
    B = B/max(B(:));
    Bsvg{i} = B;
    if isnan(max(B(:)))
        warning('NaN problem.');
    end
end
 
%%
% Do the rendering.

col = [ [1 0 0]; [0 1 0]; [0 0 1]; [1 1 0]; [1/2 1 1/2]; [1/2 1/2 1] ];

for i=1:Q
    w = W(:,i)'; w = w/sum(w);
    B = Bsvg{i};
    % display
    clf; 
    opts.alpha = 1; % transparency
	opts.color = sum(col(1:p,:) .* repmat(w(:),[1 3]) );
    opts.isolevel = median(B(:));
    opts.isolevel = (max(B(:))-min(B(:)))/2;
    F = plot_isosurface(B,opts);    
    % save as image
    if p==4
        str = ['barycenter-' num2str(S(i)*(q-1)) '-' num2str(T(i)*(q-1))];
    else
        str = ['barycenter-' num2str(i)];
    end
    saveas(gcf, [rep str '.png'], 'png');
    % save as object
    if 0
    mesh.vertices = F.vertices;
    mesh.triangles = F.faces;
    mesh.numVertices = size(F.vertices, 1);
    vertexTexture = zeros(mesh.numVertices,1);
    writeTexturedObj([rep 'barycenter-' str '.obj'], mesh, vertexTexture); 
    end
end
