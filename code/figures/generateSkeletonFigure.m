%% Load mesh

%[X,T] = readOff('../data/meshes/meshr05.off');
% [X,T] = readOff('../data/meshes/double-torus.off');
%[X,T] = readOff('../data/meshes/108_simplified.off');
[X,T] = readOff('../data/meshes/wolf0.off');
M = getMeshData(X,T,10); % compute 10 LB eigenfunctions for fun

%% Make human graph

% Human
%edges = [1 2; 2 6; 3 4; 4 5; 5 6; 6 7; 7 8; 8 9; 6 10; 10 11; 11 12; 12 13; 11 14; 14 15];

% For double torus
% edges = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 1; 2 9; 9 6];

% For chair
% edges = [1 5; 2 7; 4 9; 3 11; ... % legs
%     5 6; 6 7; 7 8; 8 9; 9 10; 10 11; 11 14; 14 12; 12 13; 13 5; ... % chair perimeter
%     %6 12; 19 12; 10 12; 8 12; ... % chair seat
%     5 15; 13 16; 14 17; 11 18; ... % chair back
%     15 16; 16 17; 17 18];

% For spiral
% edges = [(1:29)' (2:30)'];

% For wolf
edges = [1 2;2 3; 3 4; 3 5; 3 6; 6 7;7 8;7 9;7 10; 10 11];

edges = subdivideGraphEdges(edges); % for human
edges = subdivideGraphEdges(edges); % for double torus, chair, wolf

% Human
% fixedGraphVerts = [1;3;9;13;15];
% fixedMeshVerts = [3244;1306;1373;17;97];

% Double torus
% fixedGraphVerts = [2;6;8;4];
% fixedMeshVerts = [2080;2072;2051;2104];

% Chair
% fixedGraphVerts =  [1;    2;     3;   4;    15;   18];
% %fixedMeshVerts =  [6299; 10069; 768; 4761; 7938; 633]; %full chair
% fixedMeshVerts = [1041; 1422; 113; 852; 1187; 116];

% Spiral
% fixedGraphVerts = [1;30];
% fixedMeshVerts = [316;319];

% Wolf
fixedGraphVerts = [1;4;5;8;9;11];
fixedMeshVerts = [4343;3544;4102;1040;1453;2758];

%% Find interior points

blurTime = .005;
blurSteps = 3;
blur = @(x) blurOnMesh(x,M,blurTime,blurSteps); % faster than pre-factored?

n = M.numVertices;

% Guess some point in the interior by blurring mesh vertices
interiorPoints = zeros(length(fixedMeshVerts),3);

close all;
fig = showDescriptor(M,X(:,1));
alpha(.5);
hold on;

for i=1:length(fixedMeshVerts)
    ee = zeros(n,1);
    ee(fixedMeshVerts(i)) = 1/M.areaWeights(fixedMeshVerts(i));
    interiorPoints(i,:) = sum(bsxfun(@times,X,M.areaWeights.*blur(ee)),1);
    plot3(interiorPoints(i,1),interiorPoints(i,2),interiorPoints(i,3),'x');
    hold on
end

%% Calculate coordinates

coords = zeros(n,length(fixedMeshVerts));

for i=1:length(fixedMeshVerts)
    mvc = meanValueCoordinates(M,interiorPoints(i,:));
    mvc(mvc<0) = 0;
    mvc = mvc / sum(mvc);
    coords(:,i) = mvc./M.areaWeights;
end

%% Compute entropy

coords(coords < 1e-15) = 1e-15;
entropies = -sum(bsxfun(@times,coords.*log(coords),M.areaWeights));
maxEntropy = max(entropies);
targetEntropy = maxEntropy;

%% Solve propagation problem

graphEdges = [edges ; edges(:,2) edges(:,1)];
edgeWeights = ones(1,size(graphEdges,1))';
fixedVertices = fixedGraphVerts;
fixedDists = coords;
areaWeights = M.areaWeights;

blurTime = .0015; % .001 for humans, double torus
blurSteps = 3;
kernel = @(x) blurOnMesh(x,M,blurTime,blurSteps); 
kernelTranspose = @(x) blurOnMesh(x,M,blurTime,blurSteps,1); 

result = convolutionalPropagation(graphEdges, edgeWeights, fixedVertices, ...
                    fixedDists,areaWeights,kernel,kernelTranspose,targetEntropy);
                
%% Take expectations to find vertex positions

for i=1:size(result,2)
    graphVertexPositions(i,:) = sum(bsxfun(@times,result(:,i).*M.areaWeights,X),1);
end

close all;
fig = showDescriptor(M,X(:,1));
alpha(.5);
hold on;

for i=1:size(edges,1)
    v1 = graphVertexPositions(edges(i,1),:);
    v2 = graphVertexPositions(edges(i,2),:);
    plot3([v1(1) v2(1)],[v1(2) v2(2)],[v1(3) v2(3)],'x-');
end

%% Build up skeleton mesh

nv = max(edges(:));
sphereRadius = 4; % .05 for human
[Xsphere,Tsphere,sphereidx] = addSpheres(graphVertexPositions,[],1:nv,ones(1,nv)*sphereRadius);

sphereMesh.vertices = Xsphere;
sphereMesh.triangles = Tsphere;
sphereMesh.numVertices = size(Xsphere,1);
% showDescriptor(sphereMesh,sphereidx);

vf = graphVertexPositions(edges(:,2),:)-graphVertexPositions(edges(:,1),:);
vfBase = graphVertexPositions(edges(:,1),:);
radius = 4/3;%.01 for human, .1 for torus
nRadialSamples = 5;

[skeletonMesh,skeletonIndicator] = addCylinders(sphereMesh,vf,vfBase,radius,nRadialSamples);

skeletonIndicator(find(ismember(sphereidx,fixedGraphVerts))) = .5;

f = showDescriptor(M,X(:,1));
alpha(.1);hold on;
showDescriptor(skeletonMesh,skeletonIndicator,[],[],[],f);

%writeTexturedObj('human_skeleton2.obj', skeletonMesh, skeletonIndicator, 'skeleton.mtl');
writeTexturedObj('wolf_skeleton.obj', skeletonMesh, skeletonIndicator, 'skeleton.mtl');