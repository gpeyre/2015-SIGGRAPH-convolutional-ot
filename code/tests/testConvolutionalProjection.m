%% Set up paths

close all; clear; clc;
name = 'moomoo_s0';

%% Read shape

[X, T] = readOff([name, '.off']);
M = getMeshData(X, T);

%% Source distributions: left indicator and right indicator

x = M.vertices(:, 1);
left = zeros(M.numVertices, 1);
left(x < 0) = 1;
left = left / (left' * M.areaWeights);

right = zeros(M.numVertices, 1);
right(x >= 0) = 1;
right = right / (right' * M.areaWeights);

basis = [left right];

% Display left indicator.
showDescriptor(M, left, 'left indicator.', [], [0 2.5]);

% Display right indicator.
showDescriptor(M, right, 'right indicator.', [], [0 2.5]);

%% Target distribution.

target = zeros(M.numVertices, 1);
target(x >= -2) = 1;
target = target / (target' * M.areaWeights);

% Display target distribution.
showDescriptor(M, target, 'target distribution.', [], [0 2.5]);

%% Set up Gaussian blur function

blurTime = 1e-5; % if this gets too small, distances get noisy
blurSteps = 3;

% h = blurTime/blurSteps;
blur = @(x) applySymmetricKernel(x, M, blurTime, blurSteps);

%% Convolutional Projection. 
[convWeights] = convolutionalProjection(M.areaWeights, target, basis, blur);

%% L2 projection.
% min_{w} areaWeights' * ((target - basis * w).^2)
% s.t. sum(w) = 1
A = bsxfun(@times, basis, sqrt(M.areaWeights));
b = sqrt(M.areaWeights) .* target;
% min_{w} \|Aw - b\|^2
% s.t. sum(w) = 1
% Aw = w1A1 + w2A2 + ... wkAk = w1A1 + w2A2 + ... (1 - w1 - w2 ...)Ak
%    = Ak + w1(A1 - Ak) + w2(A2 - Ak) + ...
% Aw - b = (A1 - Ak)w1 + (A2 - Ak)w2 + ... + (Ak-1 - Ak)wk-1 - (b - Ak)
b = b - A(:, end);
A = bsxfun(@minus, A, A(:, end));
A = A(:, 1 : end - 1);
L2Weights = A \ b;
L2Weights = [L2Weights; 1 - sum(L2Weights)];

%% Display the results.
showDescriptor(M, basis * convWeights, 'convolutional projection.', [], [0 2.5]);

showDescriptor(M, basis * L2Weights, 'L2 projection.', [], [0 2.5]);
