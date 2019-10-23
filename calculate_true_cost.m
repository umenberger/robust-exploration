function cost = calculate_true_cost(A, B, ops)
% given:
% - true system A and B
% - control policy (K, Se)

% this function computes:
% - the true cost
%%
[Nx, ~] = size(B);
I = eye(Nx);
Q = ops.Q;
R = ops.R;
Se = ops.Se;
K = ops.K;
sigma_w = ops.sigma_w;

if length(sigma_w) > 1
    
    Sigma_w = diag(sigma_w.^2);
    
else
    
    Sigma_w = sigma_w^2*I;
    
end
    

ABK = A + B*K;

if spectralRadius(ABK) >= 1
    cost = inf;
else
    W = dlyap(ABK, Sigma_w + B*Se*B');
    cost = trace((Q + K'*R*K)*W) + trace(R*Se);
end
end

