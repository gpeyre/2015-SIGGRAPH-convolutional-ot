function result = convolutionalSoftCorrespondence(M1, M2, descriptorDiffs, ...
                    timeStep,lambda)

ne1 = M1.numEdges;
nv1 = M1.numVertices;
ne2 = M2.numEdges;
nv2 = M2.numVertices;

% Initialize marginals
vv1 = (ones(nv1,ne2));
ww1 = (ones(nv1,ne2));
d1 = (zeros(nv1,ne2));
vv2 = (ones(nv2,ne1));
ww2 = (ones(nv2,ne1));
d2 = (zeros(nv2,ne1));

result = (ones(nv2,nv1));

% For each vertex compute his incoming and outgoing edges
colors1 = graph_color(M1.edges);
colors2 = graph_color(M2.edges);
nc1 = max(colors1);
nc2 = max(colors2);

inColored1 = cell(nc1,1);
outColored1 = cell(nc1,1);
inColoredV1 = cell(nc1,1);
outColoredV1 = cell(nc1,1);

incomingEdges1 = cell(nv1,1);
outgoingEdges1 = cell(nv1,1);

inColored2 = cell(nc2,1);
outColored2 = cell(nc2,1);
inColoredV2 = cell(nc2,1);
outColoredV2 = cell(nc2,1);

incomingEdges2 = cell(nv2,1);
outgoingEdges2 = cell(nv2,1);

for i=1:ne1
    e1 = M1.edges(i,1);
    e2 = M1.edges(i,2);
    outgoingEdges1{e1} = [outgoingEdges1{e1} i];
    incomingEdges1{e2} = [incomingEdges1{e2} i];
    if (colors1(e1) > 0)
        outColored1{colors1(e1)} = [outColored1{colors1(e1)} i];
        outColoredV1{colors1(e1)} = [outColoredV1{colors1(e1)} e1];
    end
    if (colors1(e2) > 0)
        inColored1{colors1(e2)} = [inColored1{colors1(e2)} i];
        inColoredV1{colors1(e2)} = [inColoredV1{colors1(e2)} e2];
    end
end

for i=1:ne2
    e1 = M2.edges(i,1);
    e2 = M2.edges(i,2);
    outgoingEdges2{e1} = [outgoingEdges2{e1} i];
    incomingEdges2{e2} = [incomingEdges2{e2} i];
    if (colors2(e1) > 0)
        outColored2{colors2(e1)} = [outColored2{colors2(e1)} i];
        outColoredV2{colors2(e1)} = [outColoredV2{colors2(e1)} e1];
    end
    if (colors2(e2) > 0)
        inColored2{colors2(e2)} = [inColored2{colors2(e2)} i];
        inColoredV2{colors2(e2)} = [inColoredV2{colors2(e2)} e2];
    end
end

maxIn1 = max(cellfun(@length,incomingEdges1));
maxOut1 = max(cellfun(@length,outgoingEdges1));
maxIn2 = max(cellfun(@length,incomingEdges2));
maxOut2 = max(cellfun(@length,outgoingEdges2));

inCEdge1 = cell(nc1, maxIn1);
outCEdge1 = cell(nc1, maxOut1);
inCVtx1 = cell(nc1, maxIn1);
outCVtx1 = cell(nc1, maxOut1);
omega1 = zeros(nv1, 1);
for i=1:nv1
    if (colors1(i) > 0)
        for j=1:length(incomingEdges1{i})
            inCEdge1{colors1(i),j} = [inCEdge1{colors1(i),j} incomingEdges1{i}(j)];
            inCVtx1{colors1(i),j} = [inCVtx1{colors1(i),j} i];
        end
        for j=1:length(outgoingEdges1{i})
            outCEdge1{colors1(i),j} = [outCEdge1{colors1(i),j} outgoingEdges1{i}(j)];
            outCVtx1{colors1(i),j} = [outCVtx1{colors1(i),j} i];
        end
    end
    omega1(i) = sum(M1.edgeWeights(outgoingEdges1{i}))+sum(M1.edgeWeights(incomingEdges1{i}));
end

inCEdge2 = cell(nc2, maxIn2);
outCEdge2 = cell(nc2, maxOut2);
inCVtx2 = cell(nc2, maxIn2);
outCVtx2 = cell(nc2, maxOut2);
omega2 = zeros(nv2, 1);
for i=1:nv2
    if (colors2(i) > 0)
        for j=1:length(incomingEdges2{i})
            inCEdge2{colors2(i),j} = [inCEdge2{colors2(i),j} incomingEdges2{i}(j)];
            inCVtx2{colors2(i),j} = [inCVtx2{colors2(i),j} i];
        end
        for j=1:length(outgoingEdges2{i})
            outCEdge2{colors2(i),j} = [outCEdge2{colors2(i),j} outgoingEdges2{i}(j)];
            outCVtx2{colors2(i),j} = [outCVtx2{colors2(i),j} i];
        end
    end
    omega2(i) = sum(M2.edgeWeights(outgoingEdges2{i}))+sum(M2.edgeWeights(incomingEdges2{i}));
end

% descriptorDiffs has dimension (nv2, nv1) right now
vDesc = zeros(nv2,ne1);
wDesc = zeros(nv2,ne1);
for e=1:ne1
    v = M1.edges(e,1);
    w = M1.edges(e,2);
    vDesc(:,e) = exp(-lambda * descriptorDiffs(:,v) / (2 * timeStep * M1.edgeWeights(e) * (length(outgoingEdges1{v}) + length(incomingEdges1{v}))));
    wDesc(:,e) = exp(-lambda * descriptorDiffs(:,w) / (2 * timeStep * M1.edgeWeights(e) * (length(outgoingEdges1{w}) + length(incomingEdges1{w}))));
end

for j=1:30
    oldResult = result;
    
    for c1=1:nc1
        o1 = outColored1{c1};
        i1 = inColored1{c1};
        vtx1 = find(colors1 == c1);
        rev1 = full(sparse(vtx1,1,1:length(vtx1),nv1,1));
        for c2=1:nc2
            o2 = outColored2{c2};
            i2 = inColored2{c2};
            vtx2 = find(colors2 == c2);
            rev2 = full(sparse(vtx2,1,1:length(vtx2),nv2,1));
            
            VVO = gpuArray(vv2(:,o1));
            WWO = gpuArray(ww2(vtx2,o1));
            DO = WWO .* wDesc(vtx2,o1) .* (M2.kernel(vtx2,:)*bsxfun(@times, VVO .* vDesc(:,o1), M2.areaWeights)); 
            DO = max(DO, 1e-50);
            d2(vtx2,o1) = gather(DO);
            clear('VVO','WWO','DO');

            VVI = gpuArray(vv2(vtx2,i1));
            WWI = gpuArray(ww2(:,i1));
            DI = VVI .* vDesc(vtx2,i1) .* (M2.kernel(vtx2,:)*bsxfun(@times, WWI .* wDesc(:,i1), M2.areaWeights));
            DI = max(DI, 1e-50);
            d2(vtx2,i1) = gather(DI);
            clear('VVI','WWI','DI');

            VVO = gpuArray(vv1(:,o2));
            WWO = gpuArray(ww1(vtx1,o2));
            DO = WWO .* (M1.kernel(vtx1,:)*bsxfun(@times, VVO, M1.areaWeights));
            DO = max(DO, 1e-50);
            d1(vtx1,o2) = gather(DO);
            clear('VVO','WWO','DO');

            VVI = gpuArray(vv1(vtx1,i2));
            WWI = gpuArray(ww1(:,i2));
            DI = VVI .* (M1.kernel(vtx1,:)*bsxfun(@times, WWI, M1.areaWeights));
            DI = max(DI, 1e-50);
            d1(vtx1,i2) = gather(DI);
            clear('VVI','WWI','DI');
            
            pAll = ones(length(vtx2),length(vtx1));

            for idx=1:maxIn1
                if ~isempty(inCEdge1{c1,idx})
                    DI = gpuArray(d2(vtx2,inCEdge1{c1,idx}));
                    P = gpuArray(pAll(:,rev1(inCVtx1{c1,idx})));
                    P = P .* DI .^ bsxfun(@rdivide,M1.edgeWeights(inCEdge1{c1,idx})',bsxfun(@plus,omega2(vtx2),omega1(inCVtx1{c1,idx})'));
                    pAll(:, rev1(inCVtx1{c1,idx})) = gather(P);
                end
            end

            for idx=1:maxOut1
                if ~isempty(outCEdge1{c1,idx})
                    DO = gpuArray(d2(vtx2,outCEdge1{c1,idx}));
                    P = gpuArray(pAll(:,rev1(outCVtx1{c1,idx})));
                    P = P .* DO .^ bsxfun(@rdivide,M1.edgeWeights(outCEdge1{c1,idx})',bsxfun(@plus,omega2(vtx2),omega1(outCVtx1{c1,idx})'));
                    pAll(:, rev1(outCVtx1{c1,idx})) = gather(P);
                end
            end

            for idx=1:maxIn2
                if ~isempty(inCEdge2{c2,idx})
                    DI = gpuArray(d1(vtx1,inCEdge2{c2,idx}));
                    P = gpuArray(pAll(rev2(inCVtx2{c2,idx}),:)');
                    P = P .* DI .^ bsxfun(@rdivide,M2.edgeWeights(inCEdge2{c2,idx})',bsxfun(@plus,omega1(vtx1),omega2(inCVtx2{c2,idx})'));
                    pAll(rev2(inCVtx2{c2,idx}),:) = gather(P)';
                end
            end

            for idx=1:maxOut2
                if ~isempty(outCEdge2{c2,idx})
                    DO = gpuArray(d1(vtx1,outCEdge2{c2,idx}));
                    P = gpuArray(pAll(rev2(outCVtx2{c2,idx}),:)');
                    P = P .* DO .^ bsxfun(@rdivide,M2.edgeWeights(outCEdge2{c2,idx})',bsxfun(@plus,omega1(vtx1),omega2(outCVtx2{c2,idx})'));
                    pAll(rev2(outCVtx2{c2,idx}),:) = gather(P)';
                end
            end

            result(vtx2,vtx1) = pAll;
            ww2(vtx2,o1) = ww2(vtx2,o1) .* result(vtx2,outColoredV1{c1}) ./ d2(vtx2,o1);
            vv2(vtx2,i1) = vv2(vtx2,i1) .* result(vtx2,inColoredV1{c1}) ./ d2(vtx2,i1);
            ww1(vtx1,o2) = ww1(vtx1,o2) .* result(outColoredV2{c2},vtx1)' ./ d1(vtx1,o2);
            vv1(vtx1,i2) = vv1(vtx1,i2) .* result(inColoredV2{c2},vtx1)' ./ d1(vtx1,i2);
        end
    end
    
    integrals = sum(bsxfun(@times,M2.areaWeights,vv2.*vDesc.*(M2.kernel*(bsxfun(@times,ww2,M2.areaWeights).*wDesc))),1);
    vv2 = bsxfun(@rdivide,vv2,integrals);
    
    integrals = sum(bsxfun(@times,M1.areaWeights,vv1.*(M1.kernel*(bsxfun(@times,ww1,M1.areaWeights)))),1);
    vv1 = bsxfun(@rdivide,vv1,integrals);
    
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
    
    if change < 1e-8 && j > 5
        return
    end
end