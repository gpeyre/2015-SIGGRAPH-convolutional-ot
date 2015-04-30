function [distances,v,w] = convolutionalDistance(p0,p1,areaWeights,kernel,kernelTranspose,options)

% convolutionalDistance - compute entropy-OT ditance and rescaling factors
%   [distances,v,w] = convolutionalDistance(p0,p1,kernel,kernelTranspose, options)
%
% The normalized coupling reads
%       pi = diag(v)*K*diag(w)
% and it should satisfy
%       pi*a  = p1
%   and pi'*a = p0

n = size(p0,1);
options.null = 0;
niter = getoptions(options, 'niter', 100);
tol = getoptions(options, 'tol', 1e-6);
verb = getoptions(options, 'verb', 1);
displayFunction = getoptions(options, 'disp', @(x,y) disp(''));
disp_rate = getoptions(options, 'disp_rate', 10);
disp_time = getoptions(options, 'disp_time', 0);

if nargin<5 || isempty(kernelTranspose);
    kernelTranspose = kernel; % assume symmetry
end

if isempty(areaWeights)
    areaWeights = ones(n,1);
end

p0 = p0+eps;
p1 = p1+eps;

v = ones(size(p0));
w = ones(size(p1));

distances = zeros(size(v,2),1);

A = sum(areaWeights);
aw = bsxfun(@times,areaWeights,w);

for i=1:niter
    if disp_time
        tic
    end
    
    v = p1 ./ kernel(aw);
    av = bsxfun(@times,areaWeights,v);
    w = p0 ./ kernelTranspose(av);
    aw = bsxfun(@times,areaWeights,w);
    
    oldDistances = distances;
    
    ll = @(x) real(log(x));
    lv = av.*ll(v);
    lw = aw.*ll(w);
    distances = sum((log(A)*av+lv).*kernel(aw),2) + sum(av.*kernel(lw),2);
    
    change = norm(oldDistances-distances,'fro');
    
    if verb==1, fprintf('Iteration %d:  %g\n', i, change);
    elseif verb==2, progressbar(i,niter); end
    
    if isa(displayFunction,'function_handle') && mod(i,disp_rate)==1
        displayFunction(v,w); drawnow;
    end
    
    if change<tol && i > 2
        if verb==2, progressbar(niter,niter); end
        break;
    end
       
    if disp_time, toc, end
end

distances = sqrt(max(distances,0));
fprintf('\n');