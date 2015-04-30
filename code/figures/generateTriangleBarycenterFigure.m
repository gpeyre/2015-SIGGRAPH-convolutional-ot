addpath('../toolbox/');
addpath('../image_blur/');
addpath('../blur_functions//');
addpath('../convolutional_wasserstein/');
addpath('../../data/images/shapes/'); % low-res images

rep = '../results/shapes/';
if not(exist(rep))
    mkdir(rep);
end

%% Read data

if not(exist('imresize'))
    imresize = @(x,s)image_resize(x,s);    
end

targetSize = 199;

shape_list = [2 3 4];
nImages = length(shape_list);

data = zeros(targetSize*targetSize,nImages);
for i=1:nImages
    im = imread(sprintf('shape%dfilled.png', shape_list(i)));
    im = im2double(rgb2gray(im));
    
    padding = abs(size(im,2)-size(im,1));
    pad1 = floor(padding/2);
    pad2 = padding-pad1;
    
    if size(im,1) < size(im,2)
        im = [zeros(pad1,size(im,2)) ; im ; zeros(pad2,size(im,2))];
    elseif size(im,2) < size(im,1)
        im = [zeros(size(im,1),pad1) im zeros(size(im,1),pad2)];
    end
    
    im = 1-imresize(im,[targetSize targetSize]);
    im(im<0) = 0;
    
    im = im > .01;
    
    im = im + 1e-3;
    im = im / sum(im(:));
    
    % figure;imagesc(-im);axis equal;axis off;colormap gray;
    
    data(:,i) = im(:);
end

n = size(data,1);
areaWeights = ones(n,1)/n;
data = data*n;

entropies = -sum(bsxfun(@times,data.*log(data),areaWeights),1);
maxEntropy = max(entropies);

%% Set up blur

filterSize = 1.5; %was 5
imSize = [targetSize targetSize];

if exist('imfilter')
    % using image toolbox
    h = fspecial('gaussian',[1 max(imSize)],filterSize);% hsize sigma
    h = h / sum(h);
    imBlur = @(x) imfilter(imfilter(x,h,'replicate'),h','replicate');
else
    blur = load_filtering('imgaussian', targetSize);
    imBlur = @(x)blur(x,filterSize);
end

imagesc(imBlur(reshape(data(:,1),imSize)));

blurColumn = @(x) reshape(imBlur(reshape(x,imSize)),[],1);
blurAll2 = @(x) cell2mat(cellfun(blurColumn, num2cell(x,1), 'UniformOutput', false));
blurAll = @(x) blurAll2(blurAll2(x));

% blurColumn = @(x) reshape(fastBlur(reshape(x,[targetSize,targetSize]),filterSize),[],1);
% blurAll = @(x) cell2mat(cellfun(blurColumn, num2cell(x,1), 'UniformOutput', false));

imagesc(reshape(blurAll(data(:,1)),imSize));
axis equal;
axis off;

%% Compute barycenter

entropyLimit = Inf;

steps = linspace(0,1,5);% linspace(-.5,1.5,9);

averages = cell(length(steps),length(steps));
barycenters = cell(length(steps),length(steps));

options.niter = 400; %not 300
options.verb = 2;
options.tol = 1e-15;
resh = @(x)reshape(x, imSize);
options.disp = @(x)imageplot(-resh(x));

W = { ...
    [0, 0, 1] ...
    [1, 0, 3] [0, 1, 3] ...
    [1,0,1] [1,1,2] [0,1,1] ...
    [3,0,1] [2,1,1] [1,2,1] [0,3,1] ...    
    [1,0,0] [3,1,0] [1,1,0] [1,3,0] [0,1,0] ... 
    };

B = {};
for i=1:length(W)
    alpha = W{i}; alpha = alpha/sum(alpha);
    B{i} = convolutionalBarycenter(data,alpha,areaWeights,blurAll,blurAll,entropyLimit,options);
    u = rescale(-resh(B{i}));
    imwrite(u, [rep 'shape-barycenter-' num2str(i) '.png'], 'png');
    imwrite(double(u>.5), [rep 'tresh-barycenter-' num2str(i) '.png'], 'png');
end

