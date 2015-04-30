%% Read shape

[X,T] = readOff('../data/meshes/sphere_102.off');
M = getMeshData(X,T,10); % compute 10 LB eigenfunctions for fun

%% Set up Gaussian blur function

blurTime = .005; % if this gets too small, distances get noisy
blurSteps = 50;

% blur = @(x) blurOnMesh(x,M,blurTime,blurSteps)+1e-300; % faster than pre-factored?
% blurTranspose = @(x) blurOnMesh(x,M,blurTime,blurSteps,1)+1e-300;

% structure = prefactorMeshBlur(M,blurTime,blurSteps,5000);
% blur = @(x) prefactoredBlur(x,structure,0)+1e-300; 
% blurTranspose = @(x) prefactoredBlur(x,structure,1)+1e-300;

h = blurTime/blurSteps;
nv = M.numVertices;

% Sphere is small enough that we can just write down the heat kernel explicitly
blurInverse = spdiags(M.areaWeights,0,nv,nv) - h*M.cotLaplacian;
mtx = full(blurInverse) \ diag(M.areaWeights);
mtx = mtx^blurSteps;

mtx(mtx<1e-50) = 1e-50;

blur = @(x) mtx*x;
blurTranspose = @(x) mtx'*x;

%% Make boundary

[~,fixedVerts(1)] = min(X(:,1));
[~,fixedVerts(2)] = min(X(:,2));
[~,fixedVerts(3)] = min(X(:,3));
[~,fixedVerts(4)] = max(X(:,1));
[~,fixedVerts(5)] = max(X(:,2));
[~,fixedVerts(6)] = max(X(:,3));

nFixed = length(fixedVerts);

fixedDistributions = zeros(M.numVertices,nFixed);

close all;
for i=1:nFixed
    fixedDistributions(fixedVerts(i),i) = 1/M.areaWeights(fixedVerts(i));
end

fixedDistributions = fixedDistributions+1e-10;
% fixedDistributions = blur(fixedDistributions);
showDescriptor(M,fixedDistributions(:,1));

%% Compute average entropy of data

averageEntropy = -mean(sum(bsxfun(@times,M.areaWeights,(fixedDistributions.*log(fixedDistributions))),1));

%% Solve barycenter problem using this machinery

edges = [2 1; 3 1; 4 1];
edgeWeights = [1 1 1 1 1 1];

targetEntropy = averageEntropy + 1;

result = convolutionalPropagation(edges, edgeWeights, [2 3 4], ...
    fixedDistributions(:,1:3),M.areaWeights,blur,blurTranspose,targetEntropy);

barycenter = convolutionalBarycenter(fixedDistributions(:,1:3),[1 1 1],...
    M.areaWeights,blur,blurTranspose,targetEntropy);

close all

f = subplot(1,3,1);
showDescriptor(M,sum(fixedDistributions(:,1:3),2),[],[],[],f);
title('Boundary distributions');
colorbar off;

f = subplot(1,3,2);
showDescriptor(M,result(:,1),[],[],[],f);
title('Barycenter computed using graph algorithm');
colorbar off;

f = subplot(1,3,3);
showDescriptor(M,barycenter,[],[],[],f);
title('Barycenter computed using specialized barycenter algorithm')
colorbar off;

%% Solve displacement problem using this machinery

k = 9;

edges = [(1:(k-1))' (2:k)'];
edges = [edges ; edges(:,2) edges(:,1)];
edgeWeights = ones(1,size(edges,1));

targetEntropy = averageEntropy + 1.5;

result = convolutionalPropagation(edges, edgeWeights, [1 k], ...
    fixedDistributions(:,[1 6]),M.areaWeights,blur,blurTranspose,targetEntropy);

figure

for i=1:k
    f = subplot(1,k,i);
    showDescriptor(M,result(:,i),[],[],[],f);
    colorbar off;
end


%% Solve soft mapping problem

edges = [T(:,1) T(:,2); T(:,2) T(:,3); T(:,3) T(:,1)];
edges = sort(edges,2);
edges = unique(edges,'rows');
edges = [edges; edges(:,2) edges(:,1)];

edgeWeights = ones(size(edges,1),1);

targetEntropy = averageEntropy+1.5;

result = convolutionalPropagation(edges, edgeWeights, fixedVerts, ...
                    fixedDistributions,M.areaWeights,blur,blurTranspose,targetEntropy);
                
%% Show result

testVerts = randperm(M.numVertices);
nTests = 5;
testVerts = testVerts(1:nTests);

figure;
for i=1:length(testVerts)
    f = subplot(2,length(testVerts),i);
    ind = zeros(M.numVertices,1);
    ind(testVerts(i)) = 1;
    showDescriptor(M,ind,[],[],[],f);
    colorbar off;
    colormap hot;
    
    f = subplot(2,length(testVerts),i+length(testVerts));
    showDescriptor(M,result(:,testVerts(i)),[],[],[],f);
    colorbar off;
    colormap hot;
end