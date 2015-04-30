%% Read shape

[X,T] = readOff('../data/meshes/cat1.off');
M = getMeshData(X,T,10); % compute 10 LB eigenfunctions for fun

%% Design two distributions

blurTime = .0001; % if this gets too small, distances get noisy
blurSteps = 10; % was 3

frontVtx = [18962];% 15966];
backVtx = [22553];%26142 

p0 = zeros(M.numVertices,1);
for i=1:length(frontVtx)
    p0(frontVtx(i)) = .5 / M.areaWeights(frontVtx(i));
end
p0 = blurOnMesh(p0,M,blurTime,blurSteps);
showDescriptor(M,p0);

p1 = zeros(M.numVertices,1);
for i=1:length(backVtx)
    p1(backVtx(i)) = .5 / M.areaWeights(backVtx(i));
end
p1 = blurOnMesh(p1,M,blurTime,blurSteps);
showDescriptor(M,p1);

%% Set up Gaussian blur function for barycenter

blurTime = .00015; % if this gets too small, distances get noisy
blurSteps = 10;

blur = @(x) applySymmetricKernel(x,M,blurTime,blurSteps); % faster than pre-factored?
blurTranspose = @(x) blur(x);

%% Take the barycenter

p = [p0 p1];
p(p<1e-10) = 1e-10;
nFunctions = 2;

maxEntropy = max(-sum(p.*log(p).*repmat(M.areaWeights,1,2),1));

entropyLimits = [maxEntropy (maxEntropy+1) (maxEntropy+2) (maxEntropy+3) inf];
nEntropies = length(entropyLimits);

euclideanBarycenter = sum(p,2)/nFunctions;
alpha = [1 1];

options = [];
options.niter = 100; % diminishing returns after that...
%options.unit_area_projection = 1;

barycenter = zeros(M.numVertices,nEntropies);
parfor i=1:length(entropyLimits)
    barycenter(:,i) = convolutionalBarycenter(p,alpha,M.areaWeights,blur,blurTranspose,entropyLimits(i),options);
end

save entropyTest.mat

%%

for i=1:nEntropies
    showDescriptor(M,barycenter(:,i));
end

%% Write meshes

clear full
for i=1:nEntropies
    e = entropyLimits(i);
    e = full(e);
    name = sprintf('entropy_barycenter_%g.obj',e);
    writeTexturedObj(name, M, barycenter(:,i), 'scalar_function.mtl');
end

%% Write boundary

writeTexturedObj('p0.obj', M, p0, 'scalar_function.mtl');
writeTexturedObj('p1.obj', M, p1, 'scalar_function.mtl');