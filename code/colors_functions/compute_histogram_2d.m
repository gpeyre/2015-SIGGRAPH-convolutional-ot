function [H,ai,bi] = compute_histogram_2d(a,b, N, A,B)

% compute_histogram_2d - compute a 2D histogram
%
%   H = compute_histogram_2d(a,b,N; A,B);
%
%   a and b are two vectors.
%   H is an (N,N) matrices representing the histogram of values of (a,b).
%
%   A=[amin,amax], B=[bmin,bmax] indicates the range of the histograms.
%
%   Copyright (c) 2014 Gabriel Peyre

% save shapes
s = size(a); 

a = a(:);
b = b(:);

if nargin<4
    A = [min(a) max(a)];
end
if nargin<5
    B = [min(b) max(b)];
end

% rescale
a = (a-A(1))/(A(2)-A(1));
b = (b-B(1))/(B(2)-B(1));
% quantize
P = length(a);
m = 1-1e-10;
ai = clamp( 1 + floor(a*m*N), 1, N);
bi = clamp( 1 + floor(b*m*N), 1, N);

H = full( sparse(ai,bi,ones(P,1), N,N)/P );

ai = reshape(ai,s);
bi = reshape(bi,s);
