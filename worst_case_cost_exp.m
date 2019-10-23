function res = worst_case_cost_exp(ops)
% given:
% - nominal model
% - control policy
% - initial parameter uncertainty

% this function computes:
% - the worst case cost

%%
% nominal model
Ab = ops.A;
Bb = ops.B;

Kdc = ops.K;

Se = ops.Se;

[Nx,Nu] = size(Bb);

% cost function
Q = ops.Q;
R = ops.R;

% disturbance covariance
sigma_w = ops.sigma_w;

D0 = ops.D; % initial uncertainty

delta = ops.delta;

const = 1/(sigma_w*(sqrt(Nx+Nu)+sqrt(Nx)+sqrt(2*log(1/delta))))^2;  % Dean's bound

if isfield(ops,'ellipsoidal_uncertainty')
    % use bound in Proposition 2.1
    if ops.ellipsoidal_uncertainty
        
        const = 1/(sigma_w^2*chi2inv(1-delta, Nx*Nx + Nx*Nu));
        
    end
end

%% compute worst-case state covariance for given controller

yalmip('clear')

W = sdpvar(Nx,Nx);

tau = sdpvar(1); % J: changed from t and now defined as dec. variable

I = eye(Nx);
H = blkdiag(W,Se);
F = [W*Ab' + (W*Kdc')*Bb'; Se*Bb'];
G = -[W W*Kdc'; zeros(Nu,Nx) Se];

T1 = [H,  F, G;
      F', W - sigma_w^2*I, zeros(Nx,Nx+Nu);
      G', zeros(Nx+Nu,Nx), zeros(Nx+Nu)];

T2 = blkdiag(zeros(Nx+Nu,Nx+Nu), I, -const*(D0)); % D = const*(Phi'*Phi)

Objective = trace((Q+Kdc'*R*Kdc)*W) + trace(R*Se);
Constraints = [T1-tau*T2>=0, W>=0, tau>=0];
ops = sdpsettings('solver','mosek','verbose',0);
sol = optimize(Constraints,Objective, ops);

cost = double(Objective);

%%

res.sol = sol;
res.cost = cost;
res.W = double(W);
res.t = double(tau);
res.T2 = T2;
res.const = const;

res.lmi = double(T1-tau*T2);

end