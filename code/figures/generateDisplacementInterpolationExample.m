%% Read shape

[X,T] = readOff('../data/meshes/shears.off');
M = getMeshData(X,T,10); % compute 10 LB eigenfunctions for fun

MM = M;
MM.vertices(:,3) = -M.vertices(:,2);
MM.vertices(:,2) = M.vertices(:,3);
M = MM;

%% Design a few functions to average

centerVerts = [168 801 5256 5413];
nFunctions = length(centerVerts);

distributions = zeros(M.numVertices,nFunctions);

for i=1:nFunctions
    distributions(centerVerts(i),i) = 1./M.areaWeights(centerVerts(i));
end

distributions = blurOnMesh(distributions,M,.0001,3);

close all
for i=1:nFunctions
    f = subplot(1,nFunctions,i);
    showDescriptor(M,distributions(:,i),[],[],[],f);
    colorbar off;
    title(sprintf('Distribution %d',i));
end

%% Set up Gaussian blur function

blurTime = .001; % if this gets too small, distances get noisy
blurSteps = 200;

blur = @(x) applySymmetricKernel(x,M,blurTime,blurSteps); % faster than pre-factored?
blurTranspose = @(x) blur(x);

%% Try displacement interpolation

nTimeSteps = 5;

p1 = (distributions(:,1)+distributions(:,2))/2;
p2 = (distributions(:,3)+distributions(:,4))/2;

% Optional entropy bound
e1 = -sum(p1.*log(p1).*M.areaWeights);
e2 = -sum(p2.*log(p2).*M.areaWeights);
eMax = max(e1,e2); % make Inf for basic optimization
interp = zeros(M.numVertices,nTimeSteps);

for i=1:nTimeSteps
    fprintf('i = %d\n',i);
    t = (i-1)/(nTimeSteps-1);
    alpha = [t 1-t]; 
    
    options = [];
    options.niter = 100; % diminishing returns after that...
    
    options.unit_area_projection = 1;
    
    [b,v] = convolutionalBarycenter([p1 p2],alpha,M.areaWeights,blur,blurTranspose,Inf,options);
    options.initial_v = v;
    options.initial_barycenter = b;
    interp(:,i) = convolutionalBarycenter([p1 p2],alpha,M.areaWeights,blur,blurTranspose,eMax,options);
end

%% Animate the result

f = figure;
for i=1:nTimeSteps
    clf;
    showDescriptor(M,interp(:,i),[],[],[],f);
    colorbar off;
    drawnow;
    pause(0.25);
end

%% Write out images

for i=1:1:nTimeSteps
    fprintf('i = %d\n',i);
    t = (i-1)/(nTimeSteps-1);
    writeTexturedObj(sprintf('p%g_entropy_bounded.obj',t), M, interp(:,i), 'scalar_function.mtl');
end

save entropy_bounded_result.mat