%% Read big shape

[BX1,T1] = readOff('../data/meshes/meshr00.off');
M1 = getMeshData(BX1,T1,200); % compute 150 LB eigenfunctions for fun
[BX2,T2] = readOff('../data/meshes/meshr05.off');
M2 = getMeshData(BX2,T2,200); % compute 150 LB eigenfunctions for fun

nTimes = 100;
bigWks1 = waveKernelSignature(M1, nTimes);
bigWks2 = waveKernelSignature(M2, nTimes);

%% Read shape

[X1,T1] = readOff('../data/meshes/meshr00_1024.off');
M1 = getMeshData(X1,T1,10); % compute 150 LB eigenfunctions for fun
M1.areaWeights = full(M1.areaWeights);
edgeWeights = (X1(M1.edges(:,1),:) - X1(M1.edges(:,2),:)).^2;
edgeWeights = sum(edgeWeights,2);
edgeWeights = 1 ./edgeWeights;
M1.edgeWeights = edgeWeights;
[X2,T2] = readOff('../data/meshes/meshr05_1024.off');
M2 = getMeshData(X2,T2,10); % compute 150 LB eigenfunctions for fun
M2.areaWeights = full(M2.areaWeights);
edgeWeights = (X2(M2.edges(:,1),:) - X2(M2.edges(:,2),:)).^2;
edgeWeights = sum(edgeWeights,2);
edgeWeights = 1 ./edgeWeights;
M2.edgeWeights = edgeWeights;

%% Set up Gaussian blur function

blurTime = .001; % if this gets too small, distances get noisy
blurSteps = 50;

h = blurTime/blurSteps;

nv1 = M1.numVertices;

blurInverse = spdiags(M1.areaWeights,0,nv1,nv1) - h*M1.cotLaplacian;
mtx1 = full(blurInverse) \ diag(M1.areaWeights);
mtx1 = mtx1^blurSteps;
mtx1 = mtx1 * diag(1 ./ M1.areaWeights);

mtx1(mtx1<1e-50) = 1e-50;
mtx1 = gpuArray(mtx1);

M1.kernel = mtx1;

nv2 = M2.numVertices;

blurInverse = spdiags(M2.areaWeights,0,nv2,nv2) - h*M2.cotLaplacian;
mtx2 = full(blurInverse) \ diag(M2.areaWeights);
mtx2 = mtx2^blurSteps;
mtx2 = mtx2 * diag(1 ./ M2.areaWeights);

mtx2(mtx2<1e-50) = 1e-50;
mtx2 = gpuArray(mtx2);

M2.kernel = mtx2;

%% Set up descriptor differences (x value only)

wks1 = bigWks1(knnsearch(BX1,X1),:);
wks2 = bigWks2(knnsearch(BX2,X2),:);

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

%descriptorDiffs = bsxfun(@minus,X1(:,1)',X2(:,1)) .^ 2;
%descriptorDiffs = (repmat(X(:,1),1,size(X,1)) - repmat(X(:,1),1,size(X,1))') .^ 2 + (repmat(X(:,2),1,size(X,1)) - repmat(X(:,2),1,size(X,1))') .^ 2;

%% Compute average entropy of data

averageEntropy = -mean(sum(bsxfun(@times,M2.areaWeights,(fixedDists.*log(fixedDists))),1));

%% Solve soft mapping problem

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
result = convolutionalSoftCorrespondence(M1,M2,descriptorDiffs,blurTime,5000);
toc
                
%% Show result

testVerts = randperm(M1.numVertices);
nTests = 5;
testVerts = testVerts(1:nTests);
%testVerts = [testVerts 410];
%testVerts = [349 293 92];
%testVerts = fixedVertices;
sources = [ 984 ... % face
            736 ... % left bicep
            233 ... % right knee
            511 ... % midsection
            329 ... % left thigh
            92  ... % left ankle
            653 ... % right forearm
            ];
%testVerts = sources;

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
