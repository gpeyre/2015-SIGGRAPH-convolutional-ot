files = dir('../data/meshes/horses/*.off');

meshes = cell(length(files),1);
distances = cell(length(files),1);
for i=1:length(files)
    offName = files(i).name;
    
    plyName = offName;
    plyName((end-2):end) = 'ply';
    
    command = sprintf('geodesic.exe ..\\data\\meshes\\horses\\%s', plyName);
    fprintf(command);
    system(command);
    
    meshFile = sprintf('..\\data\\meshes\\horses\\%s',offName);
    
    [X,T] = readOff(meshFile);
    meshes{i} = getMeshData(X,T,10);
    
    distFile = sprintf('..\\data\\meshes\\horses\\%s_distances.txt',plyName);
    distances{i} = dlmread(distFile);
end

%% Do test

gammas = [1 1e-1 1e-2 1e-3 1e-4 1e-5];
blurSteps = 10;

convDist = cell(length(meshes),length(gammas));
rescaled = cell(length(meshes),length(gammas));

for i=2:length(meshes)
    M = meshes{i};
    [~,v] = min(distances{i});
    
    for j=1:length(gammas)
        structure = prefactorMeshBlur(M,gammas(j),blurSteps,500);
        kernel = @(x) prefactoredBlur(x,structure,0)+1e-300; 
        kernelTranspose = @(x) prefactoredBlur(x,structure,1)+1e-300;
        
        p0 = zeros(M.numVertices,M.numVertices);
        p0(v,:) = 1./M.areaWeights(v);
        p1 = spdiags(1./M.areaWeights,0,M.numVertices,M.numVertices);
        
        convDist{i,j} = convolutionalDistance(p0,p1,M.areaWeights,kernel,kernelTranspose);
        rescaled{i,j} = (convDist-min(convDist))/(max(convDist)-min(convDist));
    end
end

save validation.mat