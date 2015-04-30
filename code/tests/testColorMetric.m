%%
% Test for recovering a 2-D metric field over the a*b* plane of the L*a*b*
% domain that fits the ciede2000 perceptual metric. 

addpath('../toolbox/');
addpath('../colors_functions/');
addpath('../image_blur/');
addpath('../blur_functions//');
addpath('../convolutional_wasserstein/');
% addpath('../../data/images/colors/'); % low-res images
addpath('../../data/images/colors-big/'); % high-res images


rep = '../results/color-metric/';
if not(exist(rep))
    mkdir(rep);
end

%%
% L in [0,100], A in [-85,100], B in [-107,95]
% base fixed luminance for the metric
L = 50;

%%
% helpers

mmin = @(x)min(x(:));
mmax = @(x)max(x(:));
%
eta = .001;
DiffDir = @(ab,h)deltaE2000([L ab],[L ab]+eta*[0 h])/eta;

%%
% test

ab = [3 -80];
h = [1 1];
d = DiffDir(ab,h);

%%
% fit a single metric

M = fit_metric(@(h)DiffDir(ab,h));

%%
% Display a comparison between metric and fit.

P = 100;
tlist = linspace(0,2*pi,P+1)'; tlist(end)=[];
EuclMetric = @(M,h)sqrt( h*M*h' );

D = arrayfun(@(h1,h2)DiffDir(ab,[h1 h2]), cos(tlist), sin(tlist));
D1 = arrayfun(@(h1,h2)EuclMetric(M,[h1 h2]), cos(tlist), sin(tlist));

clf; hold on;
plot( cos(tlist).*D, sin(tlist).*D, 'b' );
plot( cos(tlist).*D1, sin(tlist).*D1, 'r:' );
axis equal;

%%
% compute some image histogram

name = 'blue-7';
n = 512;
f = rescale( load_image(name, n) );
A = colorspace(['RGB->LAB'], f);
% compute histo
N = 60;
H = compute_histogram_2d(A(:,:,2),A(:,:,3),N);
delta = .0001;
% display
clf;
imageplot(log(H+delta));


%%
% Compute a field of metrics
% A in [-85,100], B in [-107,95]

% background color image
v = 100;
a = linspace(-v,v,N);
b = linspace(-v,v,N);
[B,A] = meshgrid(b,a);
L = 80*ones(N,N);
CM = colorspace('LAB->RGB', cat(3,L,A,B));
% metric
F = fit_metric_field(DiffDir, arange,brange,N);
% inverse metric
[e1,e2,l1,l2] = perform_tensor_decomp(F);
G = perform_tensor_recomp(e1,e2,1./l1,1./l2);

options.sub = 2;
clf;
h = plot_tensor_field(F, CM, options);
saveas(gcf, [rep name '-cie-metric.png'], 'png');
clf;
h = plot_tensor_field(G, CM, options);
saveas(gcf, [rep name '-cie-inverse.png'], 'png');

