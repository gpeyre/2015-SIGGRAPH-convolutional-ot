addpath('../toolbox/');
addpath('../image_blur/');
addpath('../blur_functions//');
addpath('../convolutional_wasserstein/');
addpath('../../data/images/shapes/'); % low-res images

%% Read data

if not(exist('imresize'))
    imresize = @(x,s)image_resize(x,s);    
end

targetSize = 199;
nImages = 4;

close all
data = zeros(targetSize*targetSize,nImages);
for i=1:nImages
    im = imread(sprintf('shape%dfilled.png',i));
%    im = imread(sprintf('../data/images/smiley/smiley%d.png',i));
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
    
    figure;imagesc(-im);axis equal;axis off;colormap gray;
    
    data(:,i) = im(:);
end

n = size(data,1);
areaWeights = ones(n,1)/n;
data = data*n;

entropies = -sum(bsxfun(@times,data.*log(data),areaWeights),1);
maxEntropy = max(entropies);

%% Set up blur

filterSize = 3; %was 5
imSize = [targetSize targetSize];

close all;

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

close all;

steps = linspace(0,1,5);% linspace(-.5,1.5,9);

averages = cell(length(steps),length(steps));
barycenters = cell(length(steps),length(steps));

options.niter = 300; %not 300
options.verb = 2;
options.tol = 1e-15;

close all;
for i=1:length(steps)
    parfor j=1:length(steps) % can be parallelized
        fprintf('Test %d %d\n',i,j);
        
        s = steps(i);
        t = steps(j);
        alpha = [(1-s)*(1-t) s*(1-t) (1-s)*t s*t];
        
        averages{i,j} = sum(bsxfun(@times,data,alpha),2);
        
        barycenters{i,j} = convolutionalBarycenter(data,alpha,areaWeights,blurAll,blurAll,entropyLimit,options);
        
%         imagesc(reshape(barycenters{i,j},targetSize,[])); axis equal; axis off; colormap gray; drawnow;
    end
end

%% Visualize results

isboundary = (abs(steps) < 1e-8) | (abs(steps-1) < 1e-8);

result = [];
result2 = [];

for i=1:length(steps)
    curRow = [];
    curRow2 = [];
    for j=1:length(steps)
        im = reshape(barycenters{i,j},[targetSize,targetSize]);
        im = im / max(im(:));
        
        if isboundary(i) && isboundary(j)
            im(:,1) = 1;
            im(:,end) = 1;
            im(1,:) = 1;
            im(end,:) = 1;
        end
        
        curRow = [curRow im];
        
        im = reshape(averages{i,j},[targetSize,targetSize]);
        im = im / max(im(:));
        
        if isboundary(i) && isboundary(j)
            im(:,1) = 1;
            im(:,end) = 1;
            im(1,:) = 1;
            im(end,:) = 1;
        end
        
        curRow2 = [curRow2 im];
    end
    result = [result ; curRow];
    result2 = [result2; curRow2];
end

close all
divider = 1.8;
output1 = result<max(result(:))/divider;
imagesc(output1); axis equal; axis off; colormap gray;
figure;imagesc(-result); axis equal; axis off; colormap gray;

% output2 = result2<max(result2(:))/divider;
% figure;imagesc(output2); axis equal; axis off; colormap gray;
figure;imagesc(-result2); axis equal; axis off; colormap gray;

imwrite(output1,'wasserstein.png');
% imwrite(output2,'euclideanfaces.png');
imwrite(1-result,'wasserstein_gray.png');
imwrite(1-result2,'euclidean_gray.png');
save starexample
