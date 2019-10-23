function ops_new = update_model(A, B, ops)
% given:
% - true system
% - control policy
% - nominal model

% this function:
% - updates the nominal model
%%

ops_new = ops;
Q = ops.Q;
R = ops.R;
Phi = ops.Phi;
y = ops.y;
Tss = ops.Tss;
[Nx, Nu] = size(ops.B);
Texp = ops.Texp;
xs = zeros(Nx,Texp);
us = zeros(Nu,Texp);
c = 0;

if isfield(ops,'W')
    W = ops.W;
else
    W = ops.sigma_w*randn(Nx,Texp);
end

% control policy
Ksim = ops.K;
Sesim = ops.Se;
if min(eig(Sesim)) <= 0
    Sesim = zeros(Nu, Nu);
end
for t = 1:Texp
    us(:,t) = Ksim*xs(:,t) + mvnrnd(zeros(Nu,1),Sesim,1)';
    xs(:,t+1) = A*xs(:,t) + B*us(:,t) + W(:,t);
    c = c + xs(:,t)'*Q*xs(:,t) + us(:,t)'*R*us(:,t);
end

XU = [Phi;[xs(:,1:Tss:Texp-1)',us(:,1:Tss:Texp-1)']];
x_check = [xs(:,1:Tss:Texp-1)',us(:,1:Tss:Texp-1)'];
D1_ = XU'*XU;
Y = [y;xs(:,2:Tss:Texp)'];
theta = (XU'*XU+1e-5*eye(Nx+Nu))\XU'*Y;
Ab_ = theta(1:Nx,:)';
Bb_ = theta(Nx+1:Nx+Nu,:)';

if(ops.change_estimates)
    ops_new.A = Ab_;
    ops_new.B = Bb_;
end

ops_new.exploration = c;
ops_new.D = D1_;
ops_new.y = Y;
ops_new.Phi = XU;
ops_new.check = x_check'*x_check;
end

