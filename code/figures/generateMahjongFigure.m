%% Read data

p0im = '../data/images/mahjong/eight.png';
p1im = '../data/images/mahjong/nine.png';

p0 = 1-im2double(rgb2gray(imread(p0im)));
p1 = 1-im2double(rgb2gray(imread(p1im)));
% p0 = (p0 - min(p0(:)))/(max(p0(:))-min(p0(:)));
% p1 = (p1 - min(p1(:)))/(max(p1(:))-min(p1(:)));
p0 = double(p0 > .5);
p1 = double(p1 > .5);

imSize = size(p0);
p0 = p0(:);
p1 = p1(:);

p0 = p0 + 1e-7;
p0 = p0 / sum(p0);
p1 = p1 + 1e-7;
p1 = p1 / sum(p1);

n = length(p0);
p = [p0 p1]*n;
areaWeights = ones(n,1)/n;

%%
help imfilter

%% Set up blur

filterSize = 4; %was 5

h = fspecial('gaussian',[1 max(imSize)],filterSize);% hsize sigma
h = h / sum(h);

imBlur = @(x) imfilter(imfilter(x,h,'replicate'),h','replicate');
imagesc(imBlur(reshape(p0,imSize)));

%blurColumn = @(x) reshape(fastBlur(reshape(x,imSize),filterSize),[],1);
blurColumn = @(x) reshape(imBlur(reshape(x,imSize)),[],1);
blurAll = @(x) cell2mat(cellfun(blurColumn, num2cell(x,1), 'UniformOutput', false));

imagesc(reshape(blurAll(p1),imSize));
axis equal;
axis off;

%% Compute barycenter

steps = linspace(0,1,5);
averages = cell(length(steps),1);
barycenters = cell(length(steps),1);

entropies = -sum(p.*log(p).*[areaWeights areaWeights]);
minEntropy = min(entropies);
targetEntropy = minEntropy;

close all;
for i=1:length(steps)
    fprintf('Test %d\n',i);
    t = steps(i);
    alpha = [t (1-t)]

    a = sum(bsxfun(@times,p,alpha),2);
    averages{i} = reshape(a,imSize);
    
    subplot(1,2,1);
    imagesc(averages{i});axis equal;axis off;
    colormap gray
    
    b = convolutionalBarycenter(p,alpha,areaWeights,blurAll,blurAll,targetEntropy);
    barycenters{i} = reshape(b,imSize);
    
    subplot(1,2,2);
    imagesc(barycenters{i});axis equal;axis off;
    colormap gray
    
    drawnow;
end

%%

close all
for i=1:length(barycenters)
    scale = 1;
    while 1
        ind = barycenters{i}>max(barycenters{i}(:))/scale;
        
        if sum(ind(:))>sum(p0(:)>2e-5)*.8 || sum(ind(:))>sum(p1(:)>2e-5)*.8
            break;
        end
        scale = scale + .0025;
    end
    
    ind = -ind;
    ind = (ind - min(ind(:)))/(max(ind(:))-min(ind(:)));
    imwrite(ind,sprintf('t%g.png',steps(i)));
%     imshow(ind);

    b = -barycenters{i};
    b = (b-min(b(:)))/(max(b(:))-min(b(:)));
    imwrite(b,sprintf('b%g.png',steps(i)));
    
    b = -averages{i};
    b = (b-min(b(:)))/(max(b(:))-min(b(:)));
    imwrite(b,sprintf('a%g.png',steps(i)));
end

%%

save mahjong.mat