%% Set up paths

clear;
name = 'sharp_sphere';
setupPaths;

%% Read shape

[X,T] = readOff([name '.off']);
M = getMeshData(X,T,10); % compute 10 LB eigenfunctions for fun

%% Carry out test

% We will test per-vertex distances with different blur
% radii/regularization.  This is the hardest convergence test since it will
% be delta-to-delta -- no overlap whatsoever!

blurTimes = {1e-1,1e-2,1e-3,1e-4,1e-5};
blurSteps = 10; % Can use the same for each test, whatev.

sourceVertex = 100;

source0 = zeros(M.numVertices,M.numVertices);
source0(sourceVertex,:) = 1./M.areaWeights(sourceVertex); % from vertex 1
target0 = spdiags(1./M.areaWeights,0,M.numVertices,M.numVertices); % to all other vertices

source0 = source0+1e-5;
target0 = target0+1e-5;

options.niter = 20;

integrals = sum(bsxfun(@times,source0,M.areaWeights),1);
source0 = bsxfun(@rdivide,source0,integrals);
integrals = sum(bsxfun(@times,target0,M.areaWeights),1);
target0 = bsxfun(@rdivide,target0,integrals);

distances = cell(length(blurTimes),1);

parfor curTest = 1:length(blurTimes)
    fprintf('Test %d of %d...\n',curTest,length(blurTimes));
    
    % Should pre-factor these blurs before we test this for speed
    kernel = @(x) applySymmetricKernel(x,M,blurTimes{curTest},blurSteps)+1e-300; 

    [distances{curTest},v,w] = convolutionalDistance(source0,target0,M.areaWeights,kernel,kernel,options);
end

%% Write a sample distance mesh

smallBlur = @(x) blurOnMesh(x,M,1e-3,5,0);

d = distances{5}(:,end);
d = smallBlur(d);

writeIsolineObj(M,d,'distance_1e-5.obj',25,.05,.05);

%% Write sample "delta function"

writeTexturedObj('delta_function.obj', M, source(:,1), 'scalar_function.mtl');

%% Write data out for plotting

for i=1:length(blurTimes)
    filename = sprintf('convergence_%g.dat',blurTimes{i});
    data = [(1:nIters)' experimentResults(:,i)*sqrt(blurTimes{i})];
    dlmwrite(filename,data,' ');
end