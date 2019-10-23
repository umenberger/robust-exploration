function res = multistep_rrl(ops)
% given:
% - nominal model
% - initial parameter uncertainty
% - time horizon

% this function computes:
% - rrl control policy for given time horizon

%%

% nominal system
Ab = ops.A;
Bb = ops.B;

[Nx,Nu] = size(Bb);

% cost function
Q = ops.Q;
R = ops.R;

% disturbance covariance
sigma_w = ops.sigma_w;

% initial uncertainty
D0 = ops.D;

delta = ops.delta;

const = 1/(sigma_w*(sqrt(Nx+Nu)+sqrt(Nx)+sqrt(2*log(1/delta))))^2;  %Dean's bound

if isfield(ops,'ellipsoidal_uncertainty')
    % use bound from Proposition 2.1
    if ops.ellipsoidal_uncertainty
        
        const = 1/(sigma_w^2*chi2inv(1-delta, Nx*Nx + Nx*Nu));
        
    end
    
end

% fixed multipliers
ts = ops.multipliers;

% length of each epoch
epoch_durations = ops.epoch_durations;

% horizon
horizon = length(epoch_durations);

% subsampling time
Tss = ops.Tss;


if isfield(ops,'multiplier_scalings')
    multiplier_scalings = ops.multiplier_scalings;
else
    multiplier_scalings = 1;
end


%%

best_cost = inf;
best_scaling = nan;

for k = 1:length(multiplier_scalings)
    
    
    m_scaling = multiplier_scalings(k);
    
    
    yalmip('clear')
    
    Di = D0; % initialize the uncertainty
    
    objective = 0;
    
    constraints = [];
    
    % keep track of decision variables for extraction later
    Ws = sdpvar(Nx,Nx,horizon);
    Zs = sdpvar(Nx+Nu,Nx+Nu,horizon);
    
    for i = 1:horizon
        
        % decision variables for controller i
        Z = sdpvar(Nx+Nu,Nx+Nu);
        W = Z(1:Nx, 1:Nx);
        
        % problem set-up 0
        % -------------------------------------------------------------------------
        I = eye(Nx);
        H = I;
        F = sigma_w*I;
        G = zeros(Nx, Nx+Nu);
        C = W-[Ab Bb]*Z*[Ab'; Bb'];
        B = Z*[Ab'; Bb'];
        S0_0 = blkdiag(H, C, -Z);
        S0_0(1:Nx, Nx+1:end) = [F G];
        S0_0(Nx+1:end, 1:Nx) = [F'; G'];
        S0_0(Nx+1:2*Nx, 2*Nx+1:end) = B';
        S0_0(2*Nx+1:end, Nx+1:2*Nx) = B;
        
        S0_1 = blkdiag(zeros(Nx, Nx), eye(Nx, Nx), -const*(Di));
        
        
        % objective (for problem i)
        % -------------------------------------------------------------------------
        
        
        objective = objective + trace(blkdiag(Q, R)*Z);
        
        
        % propagate the uncertainty
        % -------------------------------------------------------------------------
        
        if i < horizon
            
            Di = Di + (epoch_durations(i)/Tss)*Z;
            
        end
        
        % add constraints
        % -------------------------------------------------------------------------
        
        if i == 1 % for the first controller, we can search for the multiplier
            
            t0 = sdpvar(1);
            
            constraints = [constraints, S0_0 - t0*S0_1 >= 0, Z>=0, t0 >= 0];
            
        else
            
            constraints = [constraints, S0_0 - m_scaling*ts(i)*S0_1 >= 0, Z>=0];
            
        end
        
        Ws(:,:,i) = W;
        Zs(:,:,i) = Z;
        
    end
    
    
    % solve
    % -------------------------------------------------------------------------
    
    ops = sdpsettings('solver','mosek','verbose',0);
    sol = optimize(constraints,objective,ops);
    
    
    % check
    % -------------------------------------------------------------------------
    
    cost = double(objective);
    
    fprintf('\t scaling = %.2f, cost = %.5e\n',m_scaling,cost)
    
    if cost < best_cost
        
        Wstmp = double(Ws);
        Zstmp = double(Zs);
        
        best_cost = cost;
        best_scaling = m_scaling;
        
        sol_tmp = sol;
        
    end
    
end

%% extract solution

Ws = Wstmp;
Zs = Zstmp;

Knew = zeros(Nu,Nx,horizon);
Snew = zeros(Nu,Nu,horizon);

for i = 1:horizon
    
    W = double(Ws(:,:,i));
    K = (double(Zs(Nx+1:end,1:Nx,i)))/W;
    
    Knew(:,:,i) = K;
    Snew(:,:,i) = double(Zs(Nx+1:end,Nx+1:end,i)) - K*W*K';
    Snew(:,:,i) = (Snew(:,:,i) + Snew(:,:,i)')/2;
    
    Snew(:,:,i) = Snew(:,:,i) - min(min(eig(Snew(:,:,i))),-1e-9)*eye(Nu);
    
end


%%


res.Ks = Knew;
res.Ss = Snew;
res.sol = sol_tmp;
res.Zs = double(Zs);
res.cost = best_cost;
res.scaling = best_scaling;


fprintf('\n')
yalmip('clear')


















end