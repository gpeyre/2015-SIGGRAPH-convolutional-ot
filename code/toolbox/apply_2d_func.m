function y = apply_2d_func(f,x)

% apply_2d_func - apply a function to 2D arrays
%
%   y = apply_2d_func(f,x);
%
%   f should take a 2D NxN image and ouput a 2D NxN image
%   x is an N^2 x Q input
%   x is an N^2 x Q input
%
%   This function is used to adapt 2D functions (e.g. filtering) for the 
%   barycenter code that assume 1D vectors.
%
%   Copyright (c) 2014 Gabriel Peyre

N = sqrt(size(x,1));
P = size(x,2);

resh = @(x)reshape(x,[N N]);
flat = @(x)x(:);

y = zeros(size(x));
for i=1:P
   y(:,i) = flat(f(resh(x(:,i))));
end

end