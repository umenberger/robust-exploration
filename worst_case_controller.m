function res = worst_case_controller(ops)
% given:
% - nominal model
% - initial parameter uncertainty

% this function computes:
% - the nominal robust controller
% - the worst case cost

%%
% nominal model
Ab = ops.A;
Bb = ops.B;

[Nx,Nu] = size(Bb);

% cost function
Q = ops.Q;
R = ops.R;

% disturbance covariance
sigma_w = ops.sigma_w;

Qcl = chol(Q);
Rcl = chol(R);

D0 = ops.D; % initial uncertainty

delta = ops.delta;

const = 1/(sigma_w*(sqrt(Nx+Nu)+sqrt(Nx)+sqrt(2*log(1/delta))))^2;

if isfield(ops,'ellipsoidal_uncertainty')
    
    if ops.ellipsoidal_uncertainty
        
        const = 1/(sigma_w^2*chi2inv(1-delta, Nx*Nx + Nx*Nu));
        
    end
    
end

%% compute a new worst-case controller

yalmip('clear')

W = sdpvar(Nx,Nx);
Z = sdpvar(Nx,Nu,'full');
X = sdpvar(Nx+Nu,Nx+Nu);

tau = sdpvar(1); % J: changed from t and now defined as dec. variable

I = eye(Nx);
H = blkdiag(W, I);
F = [W*Ab'+Z*Bb'; sigma_w*I];
G = -[W Z; zeros(Nx, Nx+Nu)];

T1 = blkdiag(H, W, zeros(Nx+Nu, Nx+Nu));

T1(1:2*Nx, 2*Nx+1:end) = [F G];
T1(2*Nx+1:end, 1:2*Nx) = [F'; G'];

T2 = blkdiag(zeros(2*Nx, 2*Nx), I, -const*(D0));

C = [Qcl*W; Rcl*Z'];
S = [X C; C' W]; %the cost
Objective = trace(X);
Constraints = [T1-tau*T2>=0, W>=0,S>=0, tau>=0];
ops = sdpsettings('solver','mosek','verbose',0);
sol = optimize(Constraints,Objective, ops);

Knew = double(Z)'/(double(W));

cost = double(Objective);


%% report the cost of this worst-case controller

if (sol.problem)
    res.cost = inf;
else
    res.cost = cost;
    res.K = Knew;
    res.t = double(tau);
    res.W = double(W);
end




end