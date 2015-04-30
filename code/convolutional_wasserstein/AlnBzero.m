function [x] = AlnBzero(A, B, x0, options)

% AlnBzero - solve A * log(B * x) = 0.
%   [x] = AlnBzero(A, B, x0, iter, tol).
%
%   A is a (K, N) matrix.
%   B is a (N, K) matrix.
%   x0 is a (K, 1) column vector. It is the initial guess of x.
%   x is a (K, 1) column vector.
%
%   Copyright (c) 2015 Tao Du

% Options.
options.null = 0;
niter = getoptions(options, 'niter', 1500);
tol = getoptions(options, 'tol', 1e-7);

newx = x0;
i = 0;
while i < niter
    oldx = newx;
    % Compute F(x) = A * log(B * x).
    F = A * log(B * oldx);
    % Compute the Jacobian.
    J = A * bsxfun(@rdivide, B, B * oldx);

    % Damped Newton's method.
    stepsize = 1.0;
    scale = 0.5;
    direction = J \ F;
    newx = oldx - stepsize * direction;
    % Make sure newx is in the domain.
    while sum(sum(B * newx <= 0)) ~= 0
        stepsize = stepsize * scale;
        newx = oldx - stepsize * direction;
    end

    % Converge?
    if max(abs(newx - oldx)) < tol
        break;
    end
    i = i + 1;
end
x = newx;
end