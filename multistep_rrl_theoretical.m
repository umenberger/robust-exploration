%% compare theoretical worst case cost for nom, rrl and greedy
clc
close all
clear variables

%% true system

A = [1.1 0.5 0;
     0.0 0.9 0.1;
     0 -0.2 0.8];

B = [1 0; 
     0 0.1;
     2 0];

 
B = fliplr(B);

Q = diag([1,1,1]);
R = diag([0.1 1]);

Qcl = chol(Q);
Rcl = chol(R);

[Nx,Nu] = size(B);

N = 500; % number of experiments
Ts = 6;  % number of rollouts
sigma_u = 1;
sigma_w = 0.5;

delta = 0.05;


%% experiment control parameters

epoch_durations = 100*ones(1,10);

Tss = 1;

num_epochs = length(epoch_durations);

horizon_rrl = min(10,num_epochs);    

use_true_as_nom = 1;

num_trials = 100;

%% optimal controller

[K_opt,S_opt] = dlqr(A,B,Q,R,zeros(Nx,Nu));
K_opt = -K_opt;

%%

worst_costs = nan(3,num_epochs,num_trials);

for trial_index = 1:num_trials
    
    fprintf('\nBeginning trial %d\n',trial_index)

        
%% initial data

    x = cell(N, 1);
    u = cell(N, 1); 

    XU = [];
    Xp = [];
    
    XU_all = [];
    Xp_all = [];

    Z = [];
    for l = 1:N
        x{l}=zeros(Nx, Ts);
        u{l}=sigma_u*randn(Nu, Ts);
        x{1} = zeros(Nx, 1);
        for t=1:Ts-1
            x{l}(:, t+1) = A*x{l}(:,t) + B*u{l}(:,t) + sigma_w*randn(Nx,1); % J: added sigma_w
        end

        XU = [XU; [x{l}(:,Ts-1)' u{l}(:,Ts-1)']];
        Xp = [Xp; x{l}(:,Ts)'];
        
        XU_all = [XU_all; [x{l}(:,1:Ts-1)' u{l}(:,1:Ts-1)']];
        Xp_all = [Xp_all; x{l}(:,2:Ts)'];        


    end

% least square estimates using all the data    
    theta = (XU_all'*XU_all)\XU_all'*Xp_all;

    Ab = theta(1:Nx,:)';
    Bb = theta(Nx+1:Nx+Nu,:)';

    Ab0 = Ab;
    Bb0 = Bb;

    const = 1/(sigma_w*(sqrt(Nx+Nu)+sqrt(Nx)+sqrt(2*log(1/delta))))^2; % J: changed formula to be same as Dean paper
        
%% initial uncertainty

    D0 = XU'*XU;

    Drrl = D0;
    Dnom = D0;
    Dgreedy = D0;
    
%% set nominal parameters

% In this script we won't update the nominal parameters
    if use_true_as_nom
         Anom = A;
         Bnom = B;
    else
        Anom = Ab;
        Bnom = Bb;
    end

%% initial rrl control policy

    ops_init.A = Anom;
    ops_init.B = Bnom;
    ops_init.Q = Q;
    ops_init.R = R;
    ops_init.D = D0;
    ops_init.delta = delta;
    ops_init.sigma_w = sigma_w;

    res_init = worst_case_controller(ops_init);
    
    ops_nom = ops_init;
    
    res_nom = res_init;
    
    ops_rrl = ops_init;
    ops_rrl.Tss = Tss;
    ops_rrl.epoch_durations = epoch_durations;
    
    ops_greedy = ops_rrl;
    
    Ktmp = res_init.K;
    ops_rrl.Ks = repmat(Ktmp,1,1,horizon_rrl);
    ops_rrl.Ss = repmat(1e-5*eye(Nu),1,1,horizon_rrl);   

    Ks_rrl = repmat(Ktmp,1,1,horizon_rrl);
    Ss_rrl = repmat(1e-5*eye(Nu),1,1,horizon_rrl);
    
%% run for some epochs    


    for epoch_index = 1:num_epochs
        
        
        fprintf('\tepoch %d\n',epoch_index)
        
        Td = epoch_durations(epoch_index);
        
% controller design
% -------------------------------------------------------------------------  

% design nominal controller
        if epoch_index > 1
            res_nom = worst_case_controller(ops_nom);
        end
        
% design rrl controller
        ops_rrl.epoch_durations = epoch_durations(epoch_index:min(epoch_index+horizon_rrl-1,num_epochs));
        
        if epoch_index + horizon_rrl-1 <= num_epochs

% rrl            
            tmp = zeros(Nu,Nx,horizon_rrl);
            tmp(:,:,1:horizon_rrl-1) = Ks_rrl(:,:,2:end);
            tmp(:,:,horizon_rrl) = Ks_rrl(:,:,end);

            ops_rrl.Ks = tmp;
            
            tmp = zeros(Nu,Nu,horizon_rrl);
            tmp(:,:,1:horizon_rrl-1) = Ss_rrl(:,:,2:end);
            tmp(:,:,horizon_rrl) = Ss_rrl(:,:,end);

            ops_rrl.Ss = tmp;        

        else
            
            tmp = Ks_rrl(:,:,2:num_epochs-epoch_index+2);
            
            ops_rrl.Ks = tmp;

            tmp = Ss_rrl(:,:,2:num_epochs-epoch_index+2);
            
            ops_rrl.Ss = tmp;  

        end

        
       ops_rrl_ = ops_rrl;
       ops_rrl_.uncert_prop_wc = 1; 
       
       cost_rrl = multistep_lqr_cost(ops_rrl_);


        if epoch_index == num_epochs
            res_rrl = worst_case_controller(ops_rrl);
            Krrl = res_rrl.K;
            Srrl = zeros(Nu);
            
        else
            ops_rrl.multipliers = cost_rrl.ts;
            ops_rrl.multiplier_scalings = [linspace(0.1,1,10)];
            res_rrl = multistep_rrl(ops_rrl);
            Krrl = res_rrl.Ks(:,:,1);
            Srrl = res_rrl.Ss(:,:,1);   
            
            Ks_rrl = res_rrl.Ks;
            Ss_rrl = res_rrl.Ss;            
        end

% compute cost for this epoch 
% -------------------------------------------------------------------------  

        worst_costs(1,epoch_index,trial_index) = res_nom.cost*Td;
       
% rrl   
        ops_rrl_cost = ops_rrl;
        ops_rrl_cost.K = Krrl;
        ops_rrl_cost.Se = Srrl;
        
        res_rrl_cost = worst_case_cost_exp(ops_rrl_cost);
        
        worst_costs(2,epoch_index,trial_index) = res_rrl_cost.cost*Td;     
        
        
% run greedy exploration 
% -------------------------------------------------------------------------
% with noise budget adjusted to achieve the same cost

    res_greedy_nom = worst_case_controller(ops_greedy);
    
    ops_greedy.K = res_greedy_nom.K;

    if res_rrl_cost.cost < res_greedy_nom.cost
        ops_greedy.Se = zeros(Nu);
    else
        ops_greedy.budget = res_rrl_cost.cost;
        res_exp = optimal_exploration_noise(ops_greedy);
    end
    
    
    ops_rrl_cost = ops_greedy;
    ops_rrl_cost.K = res_greedy_nom.K;
    ops_rrl_cost.Se = res_exp.Se;

    res_exp_cost = worst_case_cost_exp(ops_rrl_cost);

    worst_costs(3,epoch_index,trial_index) = res_exp_cost.cost*Td;     
    

% update the uncertainty 
% -------------------------------------------------------------------------  

% no exploration
    K = res_nom.K;
    W = res_nom.W;
    tmp = [eye(Nx);K];
    dDnom = tmp*W*tmp';
    
    Dnom = Dnom + (Td/Tss)*dDnom;
    
    ops_nom.D = Dnom;
        
% rrl control
    K = Krrl;
    W = res_rrl_cost.W;
    tmp = [eye(Nx);K];
    dDrrl = tmp*W*tmp' + blkdiag(zeros(Nx),Srrl);  
    
    Drrl = Drrl + (Td/Tss)*dDrrl;
    
    ops_rrl.D = Drrl;
    
% greedy
    K = res_greedy_nom.K;
    W = res_exp_cost.W;
    tmp = [eye(Nx);K];
    dDrrl = tmp*W*tmp' + blkdiag(zeros(Nx),res_exp.Se);  
    
    Dgreedy = Dgreedy + (Td/Tss)*dDrrl;
    
    ops_greedy.D = Dgreedy;      



    end


        
        
        
end



%% save
save_results = 1;

if save_results
    
    s_date = datestr(clock,'ddmmmyyyy_HHMM');
    s_name = ['th_nom_rrl_greedy' s_date '.mat'];
    save(s_name,'-v7.3')    
    
    
end
%% worst case cost

figure
all_worst_cost = sum(worst_costs, 2);
nom_worst_cost = all_worst_cost(1, 1, :);
rrl_worst_cost = all_worst_cost(2, 1, :);
greedy_worst_cost = all_worst_cost(3, 1, :);
boxplot([nom_worst_cost(:) rrl_worst_cost(:) greedy_worst_cost(:)],'labels',{'nom','rrl','greedy'})
title('The theoretical worst-case cost')
grid on




