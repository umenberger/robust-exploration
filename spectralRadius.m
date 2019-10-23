function s = spectralRadius(A)
% given:
% - matrix A

% this function computes:
% - spectral radius of A

%%
if sum(sum(isnan(A)))
    s = nan;
else
    tmp = eig(A);
    s = max(sqrt(diag(tmp*tmp')));
end