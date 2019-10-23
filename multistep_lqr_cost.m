function res = multistep_lqr_cost(ops)

% given:
% - nominal model
% - sequence of policies
% - run times
% - initial parameter uncertainty

% this function computes:
% - the worst case cost
% - assuming the uncertainty evolves according to the nominal model

%% options

epoch_durations = ops.epoch_durations;

num_epochs = length(epoch_durations);

Ks = ops.Ks; % controllers
Ss = ops.Ss; % exploration

%  nominal model
A0 = ops.A;
B0 = ops.B;

[nx,nu] = size(B0);

% disturbance covariance
sigma_w = ops.sigma_w;
I = eye(nx);

Tss = ops.Tss;

if isfield(ops,'uncert_prop_wc')
    
    uncert_prop = ops.uncert_prop_wc;
    
else
    
    uncert_prop = 0;
    
end


% stuff for analytical gradient
Ws = zeros(nx,nx,num_epochs);
ts = zeros(1,num_epochs);
T2s = zeros(3*nx+nu,3*nx+nu,num_epochs);
Ds = zeros(nx+nu,nx+nu,num_epochs);

%%

ops_cost = ops;

costs = zeros(num_epochs,1);

Di = ops.D;

for ei = 1:num_epochs
    
    Ki = Ks(:,:,ei); % feedback gain
    Si = Ss(:,:,ei); % exploration variance
    
    % calculate worst-case cost for policy with current uncertainty
    ops_cost.K = Ki; % load the current controller
    ops_cost.Se = Si; % load current uncertainty
    ops_cost.D = Di; % load current uncertainty
    
    
    % This is the (time normalized) worst-case
    % cost of running the controller forever with this level of uncertainty
    % the true system
    %     res_cost = worst_case_cost_simple(ops_cost);
    
    res_cost = worst_case_cost_exp(ops_cost);
    
    if res_cost.sol.problem
        fprintf('epoch %d\n',ei)
        res_cost.sol.problem
        res_cost.sol
    end
    
    ts(ei) = res_cost.t;
    
    Ds(:,:,ei) = Di;
    
    % scale the cost by the epoch duration
    cost_i = res_cost.cost*epoch_durations(ei);
    
    costs(ei) = cost_i;
    
    % compute the state covariance
    
    if uncert_prop % worst-case uncertainty propagation
        
        Xi = res_cost.W;
        
    else % uncertainty propagation with nominal model
        
        ABK = A0 + B0*Ki;
        
        Xi = dlyap(ABK, sigma_w^2*I + B0*Si*B0');
        
    end
    
    % update the model uncertainty
    
    tmp = [eye(nx); Ki];
    Dlim = tmp*Xi*tmp' + blkdiag(zeros(nx),Si);
    
    Di = Di + 0.1*(epoch_durations(ei)/Tss)*Dlim; % note correction for subsampling
    
    
end

%%

res.cost = sum(costs);
res.costs = costs;

res.Ws = Ws;
res.T2s = T2s;
res.ts = ts;

res.Ds = Ds;

res.const = res_cost.const;







end