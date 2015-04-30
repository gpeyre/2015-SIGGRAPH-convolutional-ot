function result = prefactoredBlur(signal,structure,transpose)

result = signal(structure.ordering,:);

if ~transpose
    for i=1:structure.steps
        Lf0 = structure.reorderedLaplacian*result;
        
        for j=1:structure.batchsize:size(result,2)
            %fprintf('\tj = %d/%d\n',j,size(result,2));
            j2 = min(j + structure.batchsize - 1,size(result,2));
            result(:,j:j2) = result(:,j:j2)+structure.R\(structure.R'\Lf0(:,j:j2));
        end
    end
else
    for i=1:structure.steps
        cached = result;
        
        for j=1:structure.batchsize:size(result,2)
            %fprintf('\tj = %d/%d\n',j,size(result,2));
            j2 = min(j + structure.batchsize - 1,size(result,2));
            result(:,j:j2) = structure.R\(structure.R'\result(:,j:j2));
        end
        result = cached+structure.reorderedLaplacian'*result;
    end
end

result(structure.ordering,:) = result;