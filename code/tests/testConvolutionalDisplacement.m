%% Read shape

[X,T] = readOff('moomoo_s0.off');
M = getMeshData(X,T,10); % compute 10 LB eigenfunctions for fun

%% Set up Gaussian blur function

blurTime = .0001; % if this gets too small, distances get noisy
blurSteps = 5;

% h = blurTime/blurSteps;

blur = @(x) blurOnMesh(x,M,blurTime,blurSteps); % faster than pre-factored?
blurTranspose = @(x) blurOnMesh(x,M,blurTime,blurSteps,1);

% blurInverse = spdiags(M.areaWeights,0,nv,nv) - h*M.cotLaplacian;
% [R,p,S] = chol(blurInverse);
% blur = @(x) prefactoredBlur(x,R,S,M.areaWeights,blurSteps);

%% Compute distance from source vertex to all targets

sourceVtx = 1;
targetVtx = 1000;

% compute distances from delta function at a single source to all targets
% (so really, we're doing n Sinkhorn iterations in parallel...)
source = zeros(M.numVertices,1);
source(sourceVtx,1) = 1/M.areaWeights(sourceVtx);
target = zeros(M.numVertices,1);
target(targetVtx,1) = 1/M.areaWeights(targetVtx);

source =  blurOnMesh(source,M,.001,5);
target =  blurOnMesh(target,M,.001,5);

k = 10;

interpolation = convolutionalDisplacement(source,target,k,M.areaWeights,blur,blurTranspose);

%% Animate result

f = figure;
nSubsteps = 5;

for i=1:(k-1)
    for j=1:nSubsteps
        clf;
        t = (j-1)/nSubsteps;
        fn = interpolation(:,i)*(1-t)+interpolation(:,i+1)*t;
        showDescriptor(M,fn,[],[],[],f);
        drawnow;
    end
end

clf;
showDescriptor(M,interpolation(:,end),[],[],[],f);
drawnow;