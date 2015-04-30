function result = convolutionalPropagation(edges, edgeWeights, fixedVertices, ...
                    fixedDists,areaWeights,kernel,kernelTranspose,entropyLimit)

ne = size(edges,1);
nv = max(edges(:));
n = size(fixedDists,1);

% Initialize marginals
vv = ones(n,ne);
ww = ones(n,ne);
d = zeros(n,ne);

result = ones(n,nv);
result(:,fixedVertices) = fixedDists;

% For convenience, a binary array
isBoundary = zeros(nv,1);
isBoundary(fixedVertices) = 1;

% For each vertex compute his incoming and outgoing edges
incomingEdges = cell(nv,1);
outgoingEdges = cell(nv,1);

for i=1:ne
    e1 = edges(i,1);
    e2 = edges(i,2);
    outgoingEdges{e1} = [outgoingEdges{e1} i];
    incomingEdges{e2} = [incomingEdges{e2} i];
end

% No reason to project a distribution whose neighbors haven't been projected
queue = zeros(nv,1);
queue(fixedVertices) = 1;

projectEntropy = 0;

for j=1:1000
    toVisit = find(queue);
    
%     if mod(j,2)==0
%         toVisit = toVisit(end:-1:1); % reverse order each time for fun
%     end
%     toVisit = toVisit(randperm(length(toVisit))); % to avoid bias
    
    % Reset queue
    queue = zeros(nv,1);
    queue(fixedVertices) = 1;
    
    oldResult = result;
    
    for k=1:length(toVisit)
        v = toVisit(k);
        
        if isBoundary(v)
            p = result(:,v);
            
            for i=1:length(outgoingEdges{v}) %(v,w)
                e = outgoingEdges{v}(i);
                ww(:,e) = p ./ kernelTranspose(vv(:,e).*areaWeights);
            end
             
            for i=1:length(incomingEdges{v}) %(w,v)
                e = incomingEdges{v}(i);
                vv(:,e) = p ./ kernel(ww(:,e).*areaWeights);
            end
        else % not boundary
            omega = sum(edgeWeights(outgoingEdges{v}))+sum(edgeWeights(incomingEdges{v}));
            p = ones(n,1);
            
            for i=1:length(outgoingEdges{v}) %(v,w)
                e = outgoingEdges{v}(i);
                d(:,e) = ww(:,e).*kernelTranspose(vv(:,e).*areaWeights);
                d(:,e) = max(d(:,e),1e-10);
                p = p .* d(:,e).^(edgeWeights(e)/omega);
            end
            
            for i=1:length(incomingEdges{v}) %(w,v)
                e = incomingEdges{v}(i);
                d(:,e) = vv(:,e).*kernel(ww(:,e).*areaWeights);
                d(:,e) = max(d(:,e),0);
                p = p .* d(:,e).^(edgeWeights(e)/omega);
            end
            
            entropy = -sum(p.*log(p).*areaWeights);
            if nargin == 8 && entropy > entropyLimit && projectEntropy
                try % just ignore projection if it fails
                    fn = @(x) full(-sum(x*areaWeights.*((p.^x).*log(p))) - entropyLimit);
                    options = optimset('Display','none','tolfun',1e-6,'tolx',1e-6);
                    a = fzero(fn,[0 10],options);
                    p = p.^a;
                end
            end
            
            result(:,v) = p;
            
            for i=1:length(outgoingEdges{v}) %(v,w)
                e = outgoingEdges{v}(i);
                ww(:,e) = ww(:,e) .* p ./ d(:,e);
            end
            
            for i=1:length(incomingEdges{v}) %(w,v)
                e = incomingEdges{v}(i);
                vv(:,e) = vv(:,e) .* p ./ d(:,e);
            end
        end
        
        queue(edges(incomingEdges{v},:)) = 1;
        queue(edges(outgoingEdges{v},:)) = 1;
        queue(v) = 1;
    end
    
%     subplot(1,2,1);
%     imagesc(log(max(result,1e-20)));
%     colorbar;
%     title(sprintf('Iteration %d, log scale',j));
%     axis off;
%     
%     subplot(1,2,2)
%     imagesc(result);
%     colorbar;
%     title(sprintf('Iteration %d, result',j));
%     axis off;
%     
%     drawnow;
    
    change = sum(sum(bsxfun(@times,areaWeights,abs(result-oldResult))))/nv;
    change = full(change);
    
    fprintf('Iteration %d:  %g\n', j, change);
    
    if ~projectEntropy && change < 1e-3 && length(toVisit) == nv % was -3
        fprintf('Starting to project entropy.\n');
        projectEntropy = 1;
    elseif change < 1e-5 && j > 5 && projectEntropy
        return
    end
end