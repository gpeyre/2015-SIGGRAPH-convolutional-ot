%% Set up paths

path(path,'mesh_functions/');
path(path,'../data');

%% Read shape

[X,T] = readOff('../data/meshes/moomoo_s0.off');
M = getMeshData(X,T,10); % compute 10 LB eigenfunctions for fun

%% Set up Gaussian blur function

blurTime = .001; % if this gets too small, distances get noisy
blurSteps = 3;

% h = blurTime/blurSteps;

blur = @(x) blurOnMesh(x,M,blurTime,blurSteps); % faster than pre-factored?
blurTranspose = @(x) blurOnMesh(x,M,blurTime,blurSteps,1);

%% Design a few functions to average

centerVerts = [300 100 600];
nFunctions = length(centerVerts);

distributions = zeros(M.numVertices,nFunctions);

for i=1:nFunctions
    distributions(centerVerts(i),i) = 1./M.areaWeights(centerVerts(i));
end

distributions = blur(distributions);

close all
for i=1:nFunctions
    f = subplot(1,nFunctions,i);
    showDescriptor(M,distributions(:,i),[],[],[],f);
    colorbar off;
    title(sprintf('Distribution %d',i));
end

%% Take the barycenter

euclideanBarycenter = sum(distributions,2)/nFunctions;

f = subplot(1,2,1);
showDescriptor(M,euclideanBarycenter,[],[],[],f);
title('Euclidean barycenter');
colorbar off;

alpha = [1 1 1];
barycenter = convolutionalBarycenter(distributions,alpha,M.areaWeights,blur,blurTranspose);

f = subplot(1,2,2);
showDescriptor(M,barycenter,[],[],[],f);
title('Wasserstein barycenter');
colorbar off;

%% Test different entropy limits

averageEntropy = -mean(sum(bsxfun(@times,distributions.*log(distributions),M.areaWeights),1));

entropyChanges = [-1.5 -1 -.5 0 .5 1 1.5];

euclideanBarycenter = sum(distributions,2)/nFunctions;

f = subplot(1,length(entropyChanges)+1,1);
showDescriptor(M,euclideanBarycenter,[],[],[],f);
title('Euclidean barycenter');
colorbar off;

for i=1:length(entropyChanges)
    targetEntropy = averageEntropy + entropyChanges(i);
    
    alpha = [1 1 1];
    barycenter = convolutionalBarycenter(distributions,alpha,M.areaWeights,blur,blurTranspose,targetEntropy);

    f = subplot(1,length(entropyChanges)+1,i+1);
    showDescriptor(M,barycenter,[],[],[],f);
    title(sprintf('entropy < average+(%g)',entropyChanges(i)));
    colorbar off;
    
    drawnow;
end
    
%% Try displacement interpolation

nTimeSteps = 100;

p1 = distributions(:,1);
p2 = distributions(:,2);

interp = zeros(M.numVertices,nTimeSteps);
for i=1:nTimeSteps
    fprintf('i = %d/%d\n',i,nTimeSteps);
    t = (i-1)/(nTimeSteps-1);
    alpha = [t 1-t]; % is this right?
    interp(:,i) = convolutionalBarycenter([p1 p2],alpha,M.areaWeights,blur,blurTranspose,averageEntropy);
end

%% Animate the result

f = figure;
for i=1:nTimeSteps
    clf;
    showDescriptor(M,interp(:,i),[],[],[],f);
    colorbar off;
    drawnow;
end