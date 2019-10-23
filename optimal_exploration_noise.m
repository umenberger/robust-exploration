function res = optimal_exploration_noise(ops)
% given:
% - nominal model
% - initial parameter uncertainty
% - user-defined budget
% - nominal controller

% this function computes:
% - optimal exploration noise Se
% - the worst case cost
%%
% nominal model
Ab = ops.A;
Bb = ops.B;

% nominal controller
K0 = ops.K;

[Nx,Nu] = size(Bb);

% cost function
Q = ops.Q;
R = ops.R;

% user-defined budget
budget = ops.budget;

% disturbance covariance
sigma_w = ops.sigma_w;

Qcl = chol(Q);
Rcl = chol(R);

D0 = ops.D; % initial uncertainty

delta = ops.delta;

const = 1/(sigma_w*(sqrt(Nx+Nu)+sqrt(Nx)+sqrt(2*log(1/delta))))^2; % Dean's bound

if isfield(ops,'ellipsoidal_uncertainty')
    % use bound in Proposition 2.1
    if ops.ellipsoidal_uncertainty
        
        const = 1/(sigma_w^2*chi2inv(1-delta, Nx*Nx + Nx*Nu));
        
    end
    
end


% T = ops.T;

%% compute worst-case state covariance for given controller

yalmip('clear')

W = sdpvar(Nx,Nx);
X = sdpvar(Nx+Nu,Nx+Nu);

tau = sdpvar(1);

Se = sdpvar(Nu,Nu,'symmetric');

I = eye(Nx);
H = blkdiag(W, Se, I);
F = [W*Ab'+W*K0'*Bb'; Se*Bb'; sigma_w*I];
G = -[W W*K0'; [zeros(Nu,Nx) Se]; zeros(Nx, Nx+Nu)];

T1 = [H, F, G; F', W, zeros(Nx,Nx+Nu); G', zeros(Nx+Nu,Nx), zeros(Nx+Nu)];

T2 = blkdiag(zeros(2*Nx+Nu, 2*Nx+Nu), I, -const*(D0)); % D = const*(Phi'*Phi)

C = [Qcl*W; Rcl*K0*W];
S = [X C; C' W]; %the cost
Objective = -trace(Se);
Constraints = [T1-tau*T2>=0,S>=0, tau >= 0, Se >=0, trace(X) + trace(R*Se) <= budget];
ops = sdpsettings('solver','mosek','verbose',0);
sol = optimize(Constraints,Objective, ops);

Se = double(Se); % variance for additive exploration noise

Wexp = double(W);

cost = double(trace(X) + trace(R*Se));


%% report the cost of this worst-case controller

res.cost = cost;
res.Se = Se;
res.sol = sol;
res.W = Wexp;


end