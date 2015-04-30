function H1 = render_lab_histogram(H,arange,brange, dispFunc)

% render_lab_histogram - create color histogram
%
%   H1 = render_lab_histogram(H,arange,brange, dispFunc)
%
%   H is a (n,n) histogram
%   H1 is a (n,n,3) color image 
%
%   Copyright (c) 2014 Gabriel Peyre

N = size(H,1);

if nargin<4
    delta = .1/N^2;
    dispFunc = @(x)log(x+delta);
end

a = linspace(arange(1),arange(2),N);
b = linspace(brange(1),brange(2),N);
[B,A] = meshgrid(b,a);
if 1
    L = 50*ones(N,N);
    CM = colorspace('LAB->RGB', cat(3,L,A,B));
    H1 = 1 - (1-CM) .* repmat(rescale(dispFunc(H)), [1 1 3]);    
else
    L = 100*rescale(-dispFunc(H));
    H1 = colorspace('LAB->RGB', cat(3,L,A,B));
end