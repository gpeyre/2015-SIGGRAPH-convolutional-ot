function F = fit_metric_field(DiffDir, arange,brange, Q)

% fit_metric_field - fit a field of tensor
%
%    F = fit_metric_field(DiffDir, arange,brange, Q);
%
%   F is a (Q,Q,2,2) field of 2x2 SDP matrices.
%   DiffDir(x,h) is the directional value of the metric at x in direction h.  
%
%   Copyright (c) 2014 Gabriel Peyre


F = zeros(Q,Q,2,2);
alist = linspace(arange(1), arange(2), Q);
blist = linspace(brange(1), brange(2), Q);
for i=1:Q
    % progressbar(i,Q);
    for j=1:Q        
        ab = [alist(i) blist(j)];
        M = fit_metric(@(h)DiffDir(ab,h));
        F(i,j,:,:) = reshape(M, [1 1 2 2]);
    end
end

end