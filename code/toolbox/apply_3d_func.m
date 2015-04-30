function y = apply_3d_func(f,x)

% apply_3d_func - apply a function to 3D arrays
%
%   y = apply_3d_func(f,x);
%
%   f should take a 2D NxNxN image and ouput a 2D NxNxN image
%   x is an N^3 x Q input
%   x is an N^3 x Q input
%
%   This function is used to adapt 3D functions (e.g. filtering) for the 
%   barycenter code that assume 1D vectors.
%
%   Copyright (c) 2014 Gabriel Peyre

N = round( (size(x,1))^(1/3) );
P = size(x,2);

resh = @(x)reshape(x,[N N N]);
flat = @(x)x(:);

y = zeros(size(x));
for i=1:P
   y(:,i) = flat(f(resh(x(:,i))));
end

end