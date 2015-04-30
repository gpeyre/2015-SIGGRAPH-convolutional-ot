function result = blurOnGraph(signal,laplacian,time,steps)

h = time/steps;
nVertices = size(laplacian,1);

timeStep = @(x) (speye(nVertices) - h*laplacian) \ x;

result = signal;
for i=1:steps
    result = timeStep(result);
    
    % help fix numerical issues
    result(result<0) = 0;
    result = result / sum(result) * sum(signal); 
end
