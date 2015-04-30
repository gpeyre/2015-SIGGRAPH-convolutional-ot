function [projectionWeights] = convolutionalProjection(areaWeights, p0, P, kernel, kernelTranspose, options)

% convolutionalProjection - solve the following optimization problem:
%   min_{PI, projectionWeights} KL(PI | kernel)
%   s.t. PI^\top * areaWeights = p0
%        PI * areaWeights = P * projectionWeights
%
%   [projectionWeights] =
%   convolutionalProjection(areaWeights, p0, P, kernel, kernelTranspose, options);
%
%   areaWeights takes into account non-uniform grids (=[] will takes
%   ones(N, 1)).
%   p0 is a distribution(N vector) to project.
%   P is a (N, K) matrix where each P(:, i) is an input density.
%   projectWeights is a K vector.
%
%   kernel takes as input an (N, :) matrix and blurs each column.
%   kernelTranspose is the adjoint filtering (=[] assumes kernel is
%   symmetric).
%
%   Copyright (c) 2015 Tao Du

if nargin < 5 || isempty(kernelTranspose)
    kernelTranspose = kernel; % assume symmetry
end
if isempty(areaWeights)
    areaWeights = ones(size(P, 1), 1);
end

% Dimensions.
[N, K] = size(P);

% Options.
options.null = 0;
niter = getoptions(options, 'niter', 1500);
tol = getoptions(options, 'tol', 1e-7);
projectionWeights = getoptions(options, 'initial_w', ones(K, 1) / K);

v = ones(N, 1);
for i = 1 : niter
    % Project into the first constraint.
    w = p0 ./ kernelTranspose(v .* areaWeights);

    % Project into the second constraint.
    b = 1 ./ (v .* kernel(w .* areaWeights));
    A = bsxfun(@times, P, areaWeights)';
    B = bsxfun(@times, P, b);

    % Solve: F(x) = A ln(B x).
    oldWeight = projectionWeights;
    projectionWeights = AlnBzero(A, B, oldWeight);

    p = P * projectionWeights;
    v = p ./ kernel(w .* areaWeights);
    if max(abs(projectionWeights - oldWeight)) < tol
        break;
    end
end

end

