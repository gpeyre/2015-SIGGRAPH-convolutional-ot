function [f1,fCeq,H1] = perform_color_transfer(K,w0,w1, fCi, L0,L1, arange,brange)

% perform_color_transfer - transfer color over Lab space
%
%   [f1,fCeq,H1] = perform_color_transfer(K,w0,w1, fCi, L0,L1, arange,brange);
%
%   Copyright (c) 2014 Gabriel Peyre

N = length(w0);

alist = linspace(arange(1),arange(2), N);
blist = linspace(brange(1),brange(2), N);
[B,A] = meshgrid(blist,alist);


%       pi = diag(w0)*K*diag(w1)
% and it should satisfy
%       pi*1  = w0 .* K(w1) = p1
%   and pi'*1 = w1 .* K(w0) = p0



% OLD -- BAD -- Transport coupling
%   pi = diag(w1)*K*diag(w0)  
% and  
%   pi*1 = w1.*K(w0) = p1 
% and 
%   pi'*1 = w0.*K(w1) = p0

p1 = w1.*K(w0);

% compute transportation maps
% diag(1/pi*1) * pi * X = (1./p1) .* w1.*K(w0 .* X)
% (A,B) -> (A0,B0) is the transport between p0 -> p1

transport = @(x,w0,w1,p1)(w1./p1) .* K(w0.*x);
A0 = transport(A, w0,w1,p1);
B0 = transport(B, w0,w1,p1);
A0(p1==0) = NaN; B0(p1==0) = NaN;

% Apply transport
ai = fCi(:,:,1);  bi = fCi(:,:,2);
u = perform_histogram_equalization( L0,L1 );
I = ai + (bi-1)*size(A0,1);
fCeq = cat(3, u, A0(I), B0(I));

H1 = compute_histogram_2d(fCeq(:,:,2), fCeq(:,:,3),N, arange, brange);
f1 = colorspace('LAB->RGB', fCeq);
