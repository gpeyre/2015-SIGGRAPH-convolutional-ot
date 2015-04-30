%%
% Test for colorspace histograms.

addpath('../toolbox/');
addpath('../colors_functions/');
addpath('../image_blur/');
addpath('../blur_functions//');
addpath('../convolutional_wasserstein/');
% addpath('../../data/images/colors/'); % low-res images
addpath('../../data/images/colors-big/'); % high-res images


% bins for histograms
N = 120;

colsp = 'XYZ';
colsp = 'LCH';  % CIE L*C*H* (CIELCH)
colsp = 'RGB';
colsp = 'LAB'; 

name = 'castle';
name = 'river';
name = 'room';
name = 'flowers-2';
name = 'blue-7';
f = rescale( load_image(name) );

fC = colorspace(['RGB->' colsp], f);

arange = range(fC(:,:,2));
brange = range(fC(:,:,3));

% L in [0,100], A in [-85,100], B in [-107,95]
myhist = @(x)hist(x(:), 100); 
mmin = @(x)min(x(:));
mmax = @(x)max(x(:));

% compute a 2D histogram
H = compute_histogram_2d(fC(:,:,2),fC(:,:,3), N, arange,brange);

H1 = render_lab_histogram(H,arange,brange);
imageplot(H1);

