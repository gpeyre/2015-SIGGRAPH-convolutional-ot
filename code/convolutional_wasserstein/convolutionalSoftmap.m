function result = convolutionalSoftmap(edges, edgeWeights, descriptorDiffs, ...
                    areaWeights,kernel,kernelTranspose,timeStep,lambda,fixedVertices,fixedDists,entropyLimit)

ne = size(edges,1);
nv = max(edges(:));
n = size(areaWeights,1);

% Initialize marginals
vv = (ones(n,ne));
ww = (ones(n,ne));
d = (zeros(n,ne));

result = (ones(n,nv));

% For each vertex compute his incoming and outgoing edges
colors = graph_color(edges);
colors(fixedVertices) = 0;
nc = max(colors);

inColored = cell(nc,1);
outColored = cell(nc,1);
inColoredV = cell(nc,1);
outColoredV = cell(nc,1);

incomingEdges = cell(nv,1);
outgoingEdges = cell(nv,1);

for i=1:ne
    e1 = edges(i,1);
    e2 = edges(i,2);
    outgoingEdges{e1} = [outgoingEdges{e1} i];
    incomingEdges{e2} = [incomingEdges{e2} i];
    if (colors(e1) > 0)
        outColored{colors(e1)} = [outColored{colors(e1)} i];
        outColoredV{colors(e1)} = [outColoredV{colors(e1)} e1];
    end
    if (colors(e2) > 0)
        inColored{colors(e2)} = [inColored{colors(e2)} i];
        inColoredV{colors(e2)} = [inColoredV{colors(e2)} e2];
    end
end

maxIn = 0;
maxOut = 0;
for i=1:nv
    maxIn = max(maxIn, length(incomingEdges{i}));
    maxOut = max(maxOut, length(outgoingEdges{i}));
end

inCEdge = cell(nc, maxIn);
outCEdge = cell(nc, maxOut);
inCVtx = cell(nc, maxIn);
outCVtx = cell(nc, maxOut);
omega = zeros(nv, 1);
for i=1:nv
    if (colors(i) > 0)
        for j=1:length(incomingEdges{i})
            inCEdge{colors(i),j} = [inCEdge{colors(i),j} incomingEdges{i}(j)];
            inCVtx{colors(i),j} = [inCVtx{colors(i),j} i];
        end
        for j=1:length(outgoingEdges{i})
            outCEdge{colors(i),j} = [outCEdge{colors(i),j} outgoingEdges{i}(j)];
            outCVtx{colors(i),j} = [outCVtx{colors(i),j} i];
        end
    end
    omega(i) = sum(edgeWeights(outgoingEdges{i}))+sum(edgeWeights(incomingEdges{i}));
end

vDesc = zeros(n,ne);
wDesc = zeros(n,ne);
for e=1:ne
    v = edges(e,1);
    w = edges(e,2);
    vDesc(:,e) = exp(-lambda * descriptorDiffs(:,v) / (2 * timeStep * edgeWeights(e) * (length(outgoingEdges{v}) + length(incomingEdges{v}))));
    wDesc(:,e) = exp(-lambda * descriptorDiffs(:,w) / (2 * timeStep * edgeWeights(e) * (length(outgoingEdges{w}) + length(incomingEdges{w}))));
end

VDO = cell(nc,1);
WDO = cell(nc,1);
VDI = cell(nc,1);
WDI = cell(nc,1);
for c=1:nc
    AO = gpuArray(repmat(areaWeights,1,length(outColored{c})));
    AI = gpuArray(repmat(areaWeights,1,length(inColored{c})));
    VDO{c} = AO .* gpuArray(vDesc(:,outColored{c}));
    WDO{c} = gpuArray(wDesc(:,outColored{c}));
    VDI{c} = gpuArray(vDesc(:,inColored{c}));
    WDI{c} = AI .* gpuArray(wDesc(:,inColored{c}));
end

EI = cell(nc,maxIn);
EO = cell(nc,maxOut);
for c=1:nc
    for j=1:maxIn
        if ~isempty(inCEdge{c,j})
            EI{c,j} = gpuArray(repmat(edgeWeights(inCEdge{c,j})' ./ omega(inCVtx{c,j})', n, 1));
        end
    end
    for j=1:maxOut
        if ~isempty(outCEdge{c,j})
            EO{c,j} = gpuArray(repmat(edgeWeights(outCEdge{c,j})' ./ omega(outCVtx{c,j})', n, 1));
        end
    end
end

for j=1:5
    oldResult = result;
    
    for c=1:nc
        VVO = gpuArray(vv(:,outColored{c}));
        WWO = gpuArray(ww(:,outColored{c}));
        DO = arrayfun(@gpu_mult, WWO, WDO{c}, kernelTranspose(VVO .* VDO{c})); 
        d(:,outColored{c}) = gather(DO);
        clear('VVO','WWO','DO');

        VVI = gpuArray(vv(:,inColored{c}));
        WWI = gpuArray(ww(:,inColored{c}));
        DI = arrayfun(@gpu_mult, VVI, VDI{c}, kernel(WWI .* WDI{c}));
        d(:,inColored{c}) = gather(DI);
        clear('VVI','WWI','DI');
        
        pAll = ones(n,nv);
        
        for idx=1:maxIn
            if ~isempty(inCEdge{c,idx})
%                pAll(:, inCVtx{c,idx}) = pAll(:, inCVtx{c,idx}) .* d(:,inCEdge{c,idx}) .^ EI{c,idx};
                DI = gpuArray(d(:,inCEdge{c,idx}));
                P = gpuArray(pAll(:,inCVtx{c,idx}));
                P = P .* DI .^ EI{c,idx};
                pAll(:, inCVtx{c,idx}) = gather(P);
            end
        end
        
        for idx=1:maxOut
            if ~isempty(outCEdge{c,idx})
%                pAll(:, outCVtx{c,idx}) = pAll(:, outCVtx{c,idx}) .* d(:,outCEdge{c,idx}) .^ EO{c,idx};
                DO = gpuArray(d(:,outCEdge{c,idx}));
                P = gpuArray(pAll(:,outCVtx{c,idx}));
                P = P .* DO .^ EO{c,idx};
                pAll(:, outCVtx{c,idx}) = gather(P);
            end
        end
        
%        pAll = pAll ./ (repmat(areaWeights' * pAll, n, 1));
        
        toVisit = find(colors == c);
        result(:,toVisit) = pAll(:, toVisit);
        if j > 50
            count = 0;
            for i=1:length(toVisit)
                v=toVisit(i);
                p = result(:,v);
                entropy = -sum(p.*log(p).*areaWeights);
                if entropy > entropyLimit + (100 - j) / 10
                    try % just ignore projection if it fails
                        fn = @(x) full(-sum(x.*areaWeights.*((p.^x).*log(p))) - entropyLimit - (100 - j) / 10);
                        options = optimset('Display','none','tolfun',1e-6,'tolx',1e-6);
                        a = fzero(fn,[0 10],options);
%                        disp(a);
                        result(:,v) = p.^a;
                        count = count + 1;
                    end
                end
            end
            if (count > 0)
                disp(count);
            end
        end
        ww(:,outColored{c}) = ww(:,outColored{c}) .* result(:, outColoredV{c}) ./ d(:, outColored{c});
        vv(:,inColored{c}) = vv(:,inColored{c}) .* result(:, inColoredV{c}) ./ d(:, inColored{c});
    end
    
    for idx=1:length(fixedVertices)
        v = fixedVertices(idx);
        result(:,v) = fixedDists(:,idx);
        p = result(:,v);

        for i=1:length(outgoingEdges{v}) %(v,w)
            e = outgoingEdges{v}(i);
            ww(:,e) = p ./ gather(wDesc(:,e).*kernelTranspose(vv(:,e).*areaWeights.*vDesc(:,e)));
        end

        for i=1:length(incomingEdges{v}) %(w,v)
            e = incomingEdges{v}(i);
            vv(:,e) = p ./ gather(vDesc(:,e).*kernel(ww(:,e).*areaWeights.*wDesc(:,e)));
        end
    end
    
    integrals = sum(bsxfun(@times,areaWeights,vv.*vDesc.*kernel(bsxfun(@times,ww,areaWeights).*wDesc)),1);
    vv = bsxfun(@rdivide,vv,integrals);
    
    subplot(1,2,1);
    imagesc(log(max(result,1e-20)));
    colorbar;
    title(sprintf('Iteration %d, log scale',j));
    axis off;
    
    subplot(1,2,2)
    imagesc(result);
    colorbar;
    title(sprintf('Iteration %d, result',j));
    axis off;
    
    drawnow;
    
    change = norm(result-oldResult,'fro');
    fprintf('Iteration %d:  %g\n', j, change);
    
    if (isnan(change))
        count = 0;
        for i=1:nv
            if nnz(isnan(result(:,i))) > 0
                count = count + 1;
                result(:,i) = oldResult(:,i);
            end
        end
        fprintf('Warning: %d vertices reverted\n', count);
    end
    
    if change < 1e-8 && j > 5
        return
    end
end