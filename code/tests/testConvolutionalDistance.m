%% Set up paths

clear;
name = 'moomoo_s0';

%% Read shape

[X,T] = readOff([name '.off']);
M = getMeshData(X,T,10); % compute 10 LB eigenfunctions for fun

%% Set up Gaussian blur function

blurTime = .00001; % if this gets too small, distances get noisy
blurSteps = 3;

% h = blurTime/blurSteps;

blur = @(x) blurOnMesh(x,M,blurTime,blurSteps); % faster than pre-factored?
blurTranspose = @(x) blurOnMesh(x,M,blurTime,blurSteps,1);

% blurInverse = spdiags(M.areaWeights,0,nv,nv) - h*M.cotLaplacian;
% [R,p,S] = chol(blurInverse);
% blur = @(x) prefactoredBlur(x,R,S,M.areaWeights,blurSteps);

%% Compute distance from source vertex to all targets

sourceVtx = 1;

% compute distances from delta function at a single source to all targets
% (so really, we're doing n Sinkhorn iterations in parallel...)
source = zeros(M.numVertices,M.numVertices);
source(sourceVtx,:) = 1./M.areaWeights(sourceVtx);
target = spdiags(1./M.areaWeights,0,M.numVertices,M.numVertices);

[dist,w0,w1] = convolutionalDistance(source,target,M.areaWeights,blur,blurTranspose);

showDescriptor(M,dist);
title('delta-to-delta distance');