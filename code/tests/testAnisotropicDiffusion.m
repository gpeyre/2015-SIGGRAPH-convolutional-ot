%%
% test for 2-D anisotropic diffusion-based kernels

addpath('../toolbox/');
addpath('../colors_functions/');
addpath('../image_blur/');
addpath('../blur_functions//');
addpath('../convolutional_wasserstein/');
% addpath('../../data/images/colors/'); % low-res images
addpath('../../data/images/colors-big/'); % high-res images

% size of the 2-D histogram domain
N = 60;

%%
% Helpers

t = linspace(0,1,N);
[Y,X] = meshgrid(t,t);
Gaussian = @(m,sigma)exp(-((X-m(1)).^2+(Y-m(2)).^2)/(2*sigma^2));
Normalize = @(x)x/sum(x(:));

%%
% Load a metric.

metric_type = 'bump';
metric_type = 'uniform';
metric_type = 'aniso';
switch metric_type
    case 'uniform'
        M = zeros(N,N,2,2);
        M(:,:,1,1) = 1; M(:,:,2,2) = 1; 
    case 'bump'        
        sigma = .15;
        h = 1 ./ ( 1 + 5*Gaussian([.5 .5],sigma) );
        M = zeros(N,N,2,2);
        M(:,:,1,1) = h; M(:,:,2,2) = h; 
    case 'aniso'
        t = linspace(-.5,.5,N);
        [Y,X] = meshgrid(t,t);
        d = sqrt(X.^2+Y.^2);
        e1 = cat(3, X./d, Y./d);
        e2 = cat(3, -Y./d, X./d);
        eta = .2; % anisotropy factor
        M = perform_tensor_recomp(e1,e2,ones(N)*eta,ones(N));
end

% differential operator
[Blur, Delta,Grad] = blurAnisotropic(M);

u = zeros(N);
u([end/4 3*end/4], [end/4 3*end/4]) = 1;
% u([end/4 3*end/4], [end/2 end/2]) = 1; 

filtIter = 4*2;
mulist = [2 10 20 60]/2;
clf;
for i=1:length(mulist)
    mu = mulist(i);
    imageplot(Blur(u,mu,filtIter),['\mu=' num2str(mu)], 2,2,i );
end

%%
% OT barycenter

% load input densities
x = {[.2 .3] [.8 .85]};
x = {[.2 .8] [.8 .8]};
sigma = .02; H = {};
for i=1:2
    H{i} = Normalize( Gaussian(x{i}, sigma) ); 
end

% Gibbs kernel
mu = .03; filtIter = 10;
mu = 1; filtIter = 30;
mu = 5; filtIter = 5;
K = @(u)Blur(u,mu,filtIter);

options.tol = 1e-9;
options.niter = 50;
options.verb = 2;
options.disp = @(x)imageplot(-reshape(x,[N N]));
options.disp_rate = 2;

Q = 6;
tlist = linspace(0,1,Q);
Blist = {};
for i=1:Q
    t = tlist(i);
    w = [t 1-t];
    clf;
%    [B,u] = wassersteinBarycenter(H,K,w, [], options);
    Hv = reshape( cell2mat(H), [N*N 2]); % as matrix
    Kv = @(x)apply_2d_func(K,x);
    [B,u] = convolutionalBarycenter(Hv,w,[],Kv,[],[], options);
    B = reshape(B, [N N]);
    Blist{end+1} = -B;
end

clf;
imageplot(Blist);
