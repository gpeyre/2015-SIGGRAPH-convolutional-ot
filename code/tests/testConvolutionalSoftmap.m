%% Read shape

[X1,T1] = readOff('../data/meshes/meshr00_1024.off');
M1 = getMeshData(X1,T1,10); % compute 150 LB eigenfunctions for fun
[X2,T2] = readOff('../data/meshes/meshr05_1024.off');
M2 = getMeshData(X2,T2,10); % compute 150 LB eigenfunctions for fun

%% Set up Gaussian blur function

blurTime = .001; % if this gets too small, distances get noisy
blurSteps = 50;

% blur = @(x) blurOnMesh(x,M,blurTime,blurSteps)+1e-300; % faster than pre-factored?
% blurTranspose = @(x) blurOnMesh(x,M,blurTime,blurSteps,1)+1e-300;

% structure = prefactorMeshBlur(M,blurTime,blurSteps,5000);
% blur = @(x) prefactoredBlur(x,structure,0)+1e-300; 
% blurTranspose = @(x) prefactoredBlur(x,structure,1)+1e-300;

h = blurTime/blurSteps;
nv = M2.numVertices;

% Sphere is small enough that we can just write down the heat kernel explicitly
blurInverse = spdiags(M2.areaWeights,0,nv,nv) - h*M2.cotLaplacian;
mtx = full(blurInverse) \ diag(M2.areaWeights);
mtx = mtx^blurSteps;
mtx = mtx * diag(1 ./ M2.areaWeights);

mtx(mtx<1e-50) = 1e-50;
mtx = gpuArray(mtx);

blur = @(x) mtx*x;
blurTranspose = @(x) mtx'*x;

%% Set up descriptor differences (x value only)

nTimes = 100;
wks1 = waveKernelSignature(M1, nTimes);
wks2 = waveKernelSignature(M2, nTimes);

%%
% for HKS
% descriptorDiffs = zeros(M2.numVertices, M1.numVertices);
% for i=1:50
%     descriptorDiffs = descriptorDiffs + (repmat(wks2(:,i),1,M1.numVertices) - repmat(wks1(:,i),1,M2.numVertices)') .^ 2;
% end
% descriptorDiffs = full(descriptorDiffs .* repmat(M1.areaWeights',M2.numVertices,1));

% for WKS
descriptorDiffs = zeros(M2.numVertices, M1.numVertices);
for i=1:nTimes
    descriptorDiffs = descriptorDiffs + abs(repmat(wks2(:,i),1,M1.numVertices) - repmat(wks1(:,i),1,M2.numVertices)') ./ abs(repmat(wks2(:,i),1,M1.numVertices) + repmat(wks1(:,i),1,M2.numVertices)');
end
disp(max(max(descriptorDiffs)));
descriptorDiffs = full(descriptorDiffs .* repmat(M1.areaWeights',M2.numVertices,1));
disp(max(max(descriptorDiffs)));

%descriptorDiffs = (repmat(X2(:,1),1,size(X2,1)) - repmat(X2(:,1),1,size(X2,1))') .^ 2;
%descriptorDiffs = full(descriptorDiffs .* repmat(M.areaWeights',nv,1));
%descriptorDiffs = (repmat(X(:,1),1,size(X,1)) - repmat(X(:,1),1,size(X,1))') .^ 2 + (repmat(X(:,2),1,size(X,1)) - repmat(X(:,2),1,size(X,1))') .^ 2;

%% Make boundary

[~,fixedVerts(1)] = min(X2(:,1));
[~,fixedVerts(2)] = min(X2(:,2));
[~,fixedVerts(3)] = min(X2(:,3));
[~,fixedVerts(4)] = max(X2(:,1));
[~,fixedVerts(5)] = max(X2(:,2));
[~,fixedVerts(6)] = max(X2(:,3));

nFixed = length(fixedVerts);

fixedDistributions = zeros(M2.numVertices,nFixed);

close all;
for i=1:nFixed
    fixedDistributions(fixedVerts(i),i) = 1/M2.areaWeights(fixedVerts(i));
end

fixedDistributions = fixedDistributions+1e-10;
% fixedDistributions = blur(fixedDistributions);
showDescriptor(M2,fixedDistributions(:,1));

%% Compute average entropy of data

averageEntropy = -mean(sum(bsxfun(@times,M2.areaWeights,(fixedDists.*log(fixedDists))),1));

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

edges = [T1(:,1) T1(:,2); T1(:,2) T1(:,3); T1(:,3) T1(:,1)];
edges = sort(edges,2);
edges = unique(edges,'rows');
% edges = [edges; edges(:,2) edges(:,1)];

edgeWeights = (X1(edges(:,1),:) - X1(edges(:,2),:)).^2;
edgeWeights = sum(edgeWeights,2);
edgeWeights = 1 ./edgeWeights;

% humans: [RH LH RF LF chest back]
fixedVertices = [410 456 24 33 808 794];
fixedDists = zeros(1024,6);
fixedDists(427,2) = 1;
fixedDists(473,1) = 1;
fixedDists(17,4) = 1;
fixedDists(14,3) = 1;
fixedDists(789,5) = 1;
fixedDists(798,6) = 1;

% pigs: [tail lh rh lf rf back snout]
% fixedVertices = [237 343 281 65 769 645 1024];
% fixedDists = zeros(1024,7);
% fixedDists(321,1) = 1;
% fixedDists(243,2) = 1;
% fixedDists(247,3) = 1;
% fixedDists(785,4) = 1;
% fixedDists(766,5) = 1;
% fixedDists(772,6) = 1;
% fixedDists(1024,7) = 1;

% glasses: [right left bridge]
% fixedVertices = [1009 984 530];
% fixedDists = zeros(1024,3);
% fixedDists(857,1) = 1;
% fixedDists(856,2) = 1;
% fixedDists(537,3) = 1;

% armadillos: [lf rf lh rh chest r_ear tail]
% 293 = snout
% 292 = l_ear
%fixedVertices = [139 584 332 355 766 92 349];
%fixedDists = zeros(1024,7);
% fixedDists(135,1) = 1;
% fixedDists(915,2) = 1;
% fixedDists(71,3) = 1;
% fixedDists(91,4) = 1;
% fixedDists(724,5) = 1;
% fixedDists(479,6) = 1;
% fixedDists(850,1) = 1;
% fixedDists(107,2) = 1;
% fixedDists(322,3) = 1;
% fixedDists(374,4) = 1;
% fixedDists(781,5) = 1;
% fixedDists(505,6) = 1;
% fixedDists(166,2) = 1;
% fixedDists(984,1) = 1;
% fixedDists(54,4) = 1;
% fixedDists(486,3) = 1;
% fixedDists(809,5) = 1;
% %fixedDists(498,6) = 1;
% fixedDists(232,6) = 1;
% fixedDists(167,7) = 1;

fixedDists = fixedDists + 1e-20;
fixedDists = gather(blur(fixedDists));
fixedDists = bsxfun(@rdivide, fixedDists, (M2.areaWeights' * fixedDists));

%fixedVertices = [];
%fixedDists = [];

%%
% wks tau for meshr00/meshr05: 100
tic
result = convolutionalSoftmap(edges, edgeWeights, descriptorDiffs, ...
                    full(M2.areaWeights),blur,blurTranspose,blurTime,10,fixedVertices,fixedDists,averageEntropy+0.5);
toc
                
%% Show result

testVerts = randperm(M1.numVertices);
nTests = 5;
testVerts = testVerts(1:nTests);
%testVerts = [testVerts 410];
testVerts = [349 293 92];
testVerts = fixedVertices;
sources = [ 984 ... % face
            736 ... % left bicep
            233 ... % right knee
            511 ... % midsection
            329 ... % left thigh
            92  ... % left ankle
            653 ... % right forearm
            ];
testVerts = sources;

figure;
for i=1:length(testVerts)
    f = subplot(3,length(testVerts),i);
    ind = zeros(M1.numVertices,1);
    ind(testVerts(i)) = 1;
    showDescriptor(M1,ind,[],[],[],f);
    colorbar off;
    colormap hot;
    
    f = subplot(3,length(testVerts),i+length(testVerts));
    showDescriptor(M2,result(:,testVerts(i)),[],[],[],f);
    colorbar off;
    colormap hot;
    
    f = subplot(3,length(testVerts),i+2*length(testVerts));
    showDescriptor(M2,(exp(-200 * descriptorDiffs(:,testVerts(i)))),[],[],[],f);
    colorbar off;
    colormap hot;
end

%% Show result video

graph = sparse(edges(:,1), edges(:,2), 1 ./ sqrt(edgeWeights), 1024, 1024);
graph = graph + graph';
[dist path] = graphshortestpath(graph,387,349);

fig = figure;
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
for i=1:length(path)
    f = subplot(1,3,1);
    ind = zeros(M1.numVertices,1);
    ind(path(i)) = 1;
    showDescriptor(M1,ind,[],[],[],f);
     view([0 0 -1]);
%     camup([1 0 0]);
    colorbar off;
    colormap hot;
    
    f = subplot(1,3,2);
    showDescriptor(M2,exp(-100 * descriptorDiffs(:,path(i))),[],[],[],f);
     view([0 0 -1]);
%     camup([1 0 0]);
    colorbar off;
    colormap hot;
    
    f = subplot(1,3,3);
    showDescriptor(M2,result(:,path(i)),[],[],[],f);
     view([0 0 -1]);
%     camup([1 0 0]);
    colorbar off;
    colormap hot;
    
%     f = subplot(2,3,4);
%     showDescriptor(M2,straight_result(:,path(i)),[],[],[],f);
%     view([0 -1 0]);
%     camup([1 0 0]);
%     colorbar off;
%     colormap hot;
%     
%     f = subplot(2,3,5);
%     showDescriptor(M2,nodesc_result(:,path(i)),[],[],[],f);
%     view([0 -1 0]);
%     camup([1 0 0]);
%     colorbar off;
%     colormap hot;
%     
%     f = subplot(2,3,6);
%     showDescriptor(M2,nofixed_result(:,path(i)),[],[],[],f);
%     view([0 -1 0]);
%     camup([1 0 0]);
%     colorbar off;
%     colormap hot;
    
    F(i) = getframe(fig);
end
