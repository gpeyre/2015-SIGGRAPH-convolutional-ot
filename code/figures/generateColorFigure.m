%%
% Test for OT color transfer.

addpath('../toolbox/');
addpath('../colors_functions/');
addpath('../image_blur/');
addpath('../blur_functions//');
addpath('../convolutional_wasserstein/');
addpath('../../data/images/colors/'); % low-res images
addpath('../../data/images/colors-big/'); % high-res images

% #bins
N = 160;
N = 140;
N = 200;

%%
% helpers

mmin = @(x)min(x(:));
mmax = @(x)max(x(:));
normalize = @(h)h/sum(h(:));
setfigname = @(name)set(gcf, 'Name', name, 'NumberTitle','off');
% plot histograms
delta = .1/N^2;
dispHist = @(x)-log(x+delta);
dispCell = @(H)cellfun(dispHist, H, 'UniformOutput', false);
Entropy = @(x)-sum(x(x>0).*log(x(x>0)));

%%
% Load images

if not(exist('names'))
    names = {'street-3', 'street-5'};
    names = {'nature-3' 'art'};
    names = {'street-4' 'street-6'};
    % 
    names = {'street-2' 'street-7'};
    names = {'street-8' 'street-1'};
    % 
    names = {'nature-7' 'street-4'};
    names = {'nature-7' 'street-2'};
    %
    names = {'yellow-2' 'purple-1'};
    %
    names = {'street-11' 'street-10'};
    %
    names = {'red-9' 'blue-7'};
end

p = length(names);

rep = ['../results/color-transport/' names{1} '-' names{2} '/'];
if not(exist(rep))
    mkdir(rep);
end

f = {}; fC = {}; 
for i=1:p
    f{i} = rescale( load_image(names{i}) );
    fC{i} = colorspace('RGB->LAB', f{i});
end

%%
% Compute range of ab components

arange = [];
brange = [];
for i=1:p
    arange(:,end+1) = range(fC{i}(:,:,2));
    brange(:,end+1) = range(fC{i}(:,:,3));
end
arange = [min(arange(1,:)); max(arange(2,:))];
brange = [min(brange(1,:)); max(brange(2,:))];

%%
% Blur kernel

metric_type = 'percep';
metric_type = 'unif';
switch metric_type
    case 'unif'     
        mu = N/55;  
        mu = N/45;
        mu = N/40;
        blur = load_filtering('imgaussian', N);
        K = @(x)blur(x,mu);
    case 'percep'
        % load the metric
        L = 50; % base lighting
        eta = .001;
        DiffDir = @(ab,h)deltaE2000([L ab],[L ab]+eta*[0 h])/eta;
        M = fit_metric_field(DiffDir, arange,brange,N);
        [e1,e2,l1,l2] = perform_tensor_decomp(M);
        M = perform_tensor_recomp(e1,e2,1./sqrt(l1),1./sqrt(l2));
        Id = zeros(N,N,2,2); Id(:,:,1,1) = 1; Id(:,:,2,2) = 1;
        M = M + 0*Id; 
        % M = Id;
        % load the diffusion operator
        [blur, Delta,Grad] = blurAnisotropic(M);
        mu = .5; filtIter = 3;
        mu = .2; filtIter = 5;
        mu = .5/40; filtIter = 60;
        mu = .5/200; filtIter = 80;
        mu = .5/20; filtIter = 40;
        K = @(u)blur(u,mu,filtIter);
    otherwise
        error('Unknown metric type.');
end
Kv = @(x)apply_2d_func(K,x);   
   
%%
% Compute histograms.

H = {}; fCi = {};
for i=1:p    
    [H{i},ai,bi] = compute_histogram_2d(fC{i}(:,:,2),fC{i}(:,:,3),N, arange, brange);
    fCi{i} = cat(3,ai,bi);
end
imageplot( dispCell(H) );

%%
% Export colored histograms

for i=1:p
    H1 = render_lab_histogram(H{i},arange,brange);
    imwrite(rescale(H1), [rep 'density-' num2str(i) '.png'], 'png');
end

%% 
% Compute transport coupling
%   pi = diag(w1)*K*diag(w0)  
% and  
%   pi*1 = w1.*K(w0) = p1 
% and 
%   pi'*1 = w1.*K(w0) = p0

options.tol = 1e-9;
options.tol = 0;
options.niter = 250;
options.verb = 2;
options.disp = @(w0,w1)imageplot( dispHist(w0.*K(w1)) );

clf; 
[distances,w0,w1] = convolutionalDistance(H{1}, H{2}, [], K,[], options);

H1 = {}; f1 = {}; fCeq = {};
[f1{1},fCeq{1},H1{1}] = perform_color_transfer(K,w0, w1, fCi{1}, fC{1}(:,:,1),fC{2}(:,:,1), arange,brange);
[f1{2},fCeq{2},H1{2}] = perform_color_transfer(K,w1, w0, fCi{2}, fC{2}(:,:,1),fC{1}(:,:,1), arange,brange);

%%
% Save image

for i=1:p
    imwrite(rescale(f{i}), [rep names{i} '-original.jpg'], 'jpg');
%    imwrite(rescale(f1{i}), [rep names{i} '-' metric_type '-equalized.jpg'], 'jpg');
end

% display a comparison of histograms
figure(1); setfigname('Histograms');
clf;
imageplot(dispCell({H{1} H1{2}, H{2} H1{1}}), '', 2,2);
saveas(gcf, [rep names{1} '-' names{2} '-histograms.jpg'], 'jpg');

figure(2); setfigname('Images');
clf;
for i=1:2
    imageplot(f{i}, 'Original', 2,2,1+2*(i-1));
    imageplot(f1{i}, 'Equalized', 2,2,2+2*(i-1));
end
saveas(gcf, [rep names{1} '-' names{2} '-' metric_type '-images.jpg'], 'jpg');

%%
% Compute displacemet interpolation

options.disp = @(x)imageplot(dispHist(reshape(x, [N N])));

Q = 9;
tlist = linspace(0,1,Q);
% Q = 1; tlist = [1/2];

sharpening = 'none';
sharpening = 'min';
sharpening = 'arith';
sharpening = 'geom';

for i=1:Q
    t = tlist(i);
    w = [t 1-t];
    % select entropy bound
    switch sharpening
        case 'none'
            entropyLimit = [];
        case 'min'
            entropyLimit = min(Entropy(H{1}),Entropy(H{2}));
        case 'arith'
            entropyLimit = t*Entropy(H{1}) + (1-t)*Entropy(H{2});
        case 'geom'
            entropyLimit = exp( t*log(Entropy(H{1})) + (1-t)*log(Entropy(H{2})) );
        otherwise
            error('Unknown sharpening type');
    end
    % do the computation
    clf;   
    Hv = reshape( cell2mat(H), [N*N p]); % as matrix
    [B,u] = convolutionalBarycenter( Hv, w, [], Kv, [],entropyLimit, options);
    B = reshape(B, [N N]);
    u = mat2cell( reshape(u, [N N p]), N, N, ones(p,1) );    
    % save histogram
    H1 = render_lab_histogram(B,arange,brange);
    imwrite(rescale(H1), [rep 'barycenter-' metric_type '-' num2str(i) '.png'], 'png');
    % pi = diag(v{i})*K*diag(u{i})
    %   pi'*1 = u{i}.*K(v{i}) = bary
    %   pi*1 = v{i}.*K(u{i}) = H{i} ==> v{i} = H{i}./K(u{j})
    v = { H{1}./K(u{1}), H{2}./K(u{2}) };
    % compute average histogram 
    opt.histinterp = 1-t;
    L = perform_histogram_equalization( fC{1}(:,:,1),fC{2}(:,:,1), opt );
    % perform both equalization
    H1 = {}; f1 = {}; fCeq = {};
    [f1{1},fCeq{1},H1{1}] = perform_color_transfer(K, u{1}, v{1}, fCi{1}, fC{1}(:,:,1),L, arange,brange);
    [f1{2},fCeq{2},H1{2}] = perform_color_transfer(K, u{2}, v{2}, fCi{2}, fC{2}(:,:,1),L, arange,brange);
    % write
    imwrite(rescale(f1{1}), [rep names{1} '-barycenter-' metric_type '-' num2str(i) '.jpg'], 'jpg');
    imwrite(rescale(f1{2}), [rep names{2} '-barycenter-' metric_type '-' num2str(i) '.jpg'], 'jpg');
end
