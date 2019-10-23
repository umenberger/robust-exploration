%% compare empirical worst case cost and true system cost for nom, rrl and greedy
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

epoch_durations = 100*ones(1, 10);

Tss = 1;

num_epochs = length(epoch_durations);

horizon_rrl = min(10,num_epochs);    

change_nominal_estimates = 1;

num_trials = 100;

ellipsoidal_uncertainty = 1;

%% optimal controller

[K_opt,S_opt] = dlqr(A,B,Q,R,zeros(Nx,Nu));
K_opt = -K_opt;

%%
worst_costs = nan(3,num_epochs,num_trials);

epoch_costs = nan(3,num_epochs,num_trials);

uncertainty_measure = nan(3, num_epochs, num_trials);

execution_time = nan(5, num_epochs, num_trials);
% 1- worst controller
% 2 - rrl
% 3 - optimal noise
% 4 - multistep_lqr_cost
% 5 - worst-case controller for greedy
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
            x{l}(:, t+1) = A*x{l}(:,t) + B*u{l}(:,t) + sigma_w*randn(Nx,1); 
        end

        XU = [XU; [x{l}(:,Ts-1)' u{l}(:,Ts-1)']];
        Xp = [Xp; x{l}(:,Ts)'];
        
        XU_all = [XU_all; [x{l}(:,1:Ts-1)' u{l}(:,1:Ts-1)']];
        Xp_all = [Xp_all; x{l}(:,2:Ts)'];        


    end

% least square estimates using all the data    
%     theta = (XU_all'*XU_all)\XU_all'*Xp_all;
    theta = (XU'*XU)\XU'*Xp;

    Ab = theta(1:Nx,:)';
    Bb = theta(Nx+1:Nx+Nu,:)';

    Ab0 = Ab;
    Bb0 = Bb;

    const = 1/(sigma_w*(sqrt(Nx+Nu)+sqrt(Nx)+sqrt(2*log(1/delta))))^2; % J: changed formula to be same as Dean paper
        
    Phi = XU;
    y = Xp;
    
%% initial uncertainty

    D0 = XU'*XU;

    Drrl = D0;
    Dnom = D0;
   

%% initial rrl control policy

    ops_init.A = Ab;
    ops_init.B = Bb;
    ops_init.Q = Q;
    ops_init.R = R;
    ops_init.D = D0;
    ops_init.delta = delta;
    ops_init.sigma_w = sigma_w;
    ops_init.ellipsoidal_uncertainty = ellipsoidal_uncertainty;
    
    ops_init.Phi = Phi;
    ops_init.y = y;

    ops_init.Tss = Tss;
    ops_init.epoch_durations = epoch_durations;
    
    ops_init.exploration = 0;
    ops_init.change_estimates = change_nominal_estimates;    
    
    tic;
    res_init = worst_case_controller(ops_init);
    execution_time(1, 1, trial_index) = toc;
    
    if(res_init.cost>1e6)
    % if uncertainty is too high
        continue;
    end
    
    ops_nom = ops_init;
    ops_nom.Se = zeros(Nu);
    
    res_nom = res_init;
    
    ops_rrl = ops_init;
    ops_greedy = ops_init;

    Ktmp = res_init.K;
    
    Ks_rrl = repmat(Ktmp,1,1,horizon_rrl);
    Ss_rrl = repmat(1e-5*eye(Nu),1,1,horizon_rrl);
    
    Ks_greedy = repmat(Ktmp,1,1,horizon_rrl);
    Ss_greedy = repmat(1e-5*eye(Nu),1,1,horizon_rrl);    
    
%% run for some epochs    


    for epoch_index = 1:num_epochs
        
        
        fprintf('\tepoch %d\n',epoch_index)
        
        Td = epoch_durations(epoch_index);
        
% controller design
% -------------------------------------------------------------------------  

% design nominal controller
        if epoch_index > 1   
            tic;
            res_nom = worst_case_controller(ops_nom);
            execution_time(1, epoch_index, trial_index) = toc;
        end
        
% design rrl controller
        tmp = epoch_durations(epoch_index:min(epoch_index+horizon_rrl-1,num_epochs));
        ops_rrl.epoch_durations = tmp;
        ops_rrl.epoch_durations = tmp;
        ops_greedy.epoch_durations = tmp;
        
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
            
            % greedy      
            tmp = zeros(Nu,Nx,horizon_rrl);
            tmp(:,:,1:horizon_rrl-1) = Ks_greedy(:,:,2:end);
            tmp(:,:,horizon_rrl) = Ks_greedy(:,:,end);

            ops_greedy.Ks = tmp;

            tmp = zeros(Nu,Nu,horizon_rrl);
            tmp(:,:,1:horizon_rrl-1) = Ss_greedy(:,:,2:end);
            tmp(:,:,horizon_rrl) = Ss_greedy(:,:,end);

            ops_greedy.Ss = tmp; 

        else
            
            % rrl           
            tmp = Ks_rrl(:,:,2:num_epochs-epoch_index+2);
            
            ops_rrl.Ks = tmp;
            
            tmp = Ss_rrl(:,:,2:num_epochs-epoch_index+2);
            
            ops_rrl.Ss = tmp; 
            
            % greedy   
            tmp = Ks_greedy(:,:,2:num_epochs-epoch_index+2);
            
            ops_greedy.Ks = tmp;
            
            tmp = Ss_greedy(:,:,2:num_epochs-epoch_index+2);
            
            ops_greedy.Ss = tmp;             
            
        end

% Update multipliers        
% ------------------------------------------------------------------------- 

% rrl 
       ops_rrl_ = ops_rrl;
       ops_rrl_.uncert_prop_wc = 1; 
       tic;
       cost_rrl = multistep_lqr_cost(ops_rrl_);   
       execution_time(4, epoch_index, trial_index) = toc;
       
       ops_rrl.multipliers = cost_rrl.ts;
       ops_rrl.multiplier_scalings = [linspace(0.1,2,10)];

% greedy
       
       ops_greedy.multipliers = cost_rrl.ts;
       
       ops_greedy.multiplier_scalings = [linspace(0.1,1,10) linspace(1.01,1.5,5)];
       
% synthesise controllers       
% -------------------------------------------------------------------------  

        if epoch_index == num_epochs
            tic;
            res_rrl = worst_case_controller(ops_rrl);
            execution_time(2, epoch_index, trial_index)= toc;
            Krrl = res_rrl.K;
            Srrl = zeros(Nu); 
            
        else
            tic;
            res_rrl = multistep_rrl(ops_rrl);
            execution_time(2, epoch_index, trial_index)= toc;
            Krrl = res_rrl.Ks(:,:,1);
            Srrl = res_rrl.Ss(:,:,1);    
            
            Ks_rrl = res_rrl.Ks;
            Ss_rrl = res_rrl.Ss;
            
        end                       
        
    
% compute cost for this epoch, based on the true system
% -------------------------------------------------------------------------  

        ops_nom.K = res_nom.K;
       
        ops_rrl.K = Krrl;
        ops_rrl.Se = Srrl;

        true_nom = calculate_true_cost(A, B, ops_nom);
        true_rrl = calculate_true_cost(A, B, ops_rrl);

        epoch_costs(1,epoch_index,trial_index) = true_nom*Td;
        epoch_costs(2,epoch_index,trial_index) = true_rrl*Td;

% compute cost for this epoch, based on the worst-case empirical system
% -------------------------------------------------------------------------  

        ops_nom.K = res_nom.K;

        ops_rrl.K = Krrl;
        ops_rrl.Se = Srrl;
        
        % run greedy exploration 
% -------------------------------------------------------------------------
% with noise budget adjusted to achieve the same cost
    tic;
    res_exp_nom = worst_case_controller(ops_greedy);
    execution_time(5, epoch_index, trial_index)= toc;
    
    ops_greedy.K = res_exp_nom.K;
    res_rrl_cost = worst_case_cost_exp(ops_rrl);
    if res_rrl_cost.cost < res_exp_nom.cost 
        tic;
        ops_greedy.Se = zeros(Nu);
        execution_time(3, epoch_index, trial_index)= toc;
    else
        ops_greedy.budget = res_rrl_cost.cost;
        tic;
        res_exp = optimal_exploration_noise(ops_greedy);
        execution_time(3, epoch_index, trial_index)= toc;
        ops_greedy.K = res_exp_nom.K;
        ops_greedy.Se = res_exp.Se;
    end
    
    true_greedy = calculate_true_cost(A, B, ops_greedy);
    epoch_costs(3,epoch_index,trial_index) = true_greedy*Td;
    res_greedy_cost = worst_case_cost_exp(ops_greedy);

    
    worst_costs(1,epoch_index,trial_index) = res_nom.cost*Td;
    worst_costs(2,epoch_index,trial_index) = res_rrl_cost.cost*Td;
    worst_costs(3,epoch_index,trial_index) = res_greedy_cost.cost*Td;     
 
% update the uncertainty 
% -------------------------------------------------------------------------  

        W = sigma_w*randn(Nx,Td);
        
        ops_nom.W = W;
        ops_rrl.W = W;
        ops_greedy.W = W;
        
        ops_nom.Texp = Td;
        ops_rrl.Texp = Td;
        ops_greedy.Texp = Td;

        % update models
        ops_nom = update_model(A, B, ops_nom);
        ops_rrl = update_model(A, B, ops_rrl); 
        ops_greedy = update_model(A, B, ops_greedy);

        uncertainty_measure(1,epoch_index,trial_index) = trace(inv(ops_nom.D));
        uncertainty_measure(2,epoch_index,trial_index) = trace(inv(ops_rrl.D));
        uncertainty_measure(3,epoch_index,trial_index) = trace(inv(ops_greedy.D));

    end


        
        
        
end


%%    
    
fprintf('cost on true\n')
epoch_costs
sum(epoch_costs,2)

fprintf('worst case cost\n')
worst_costs
sum(worst_costs,2)

%% save

save_results = 1;

if save_results
    
    s_date = datestr(clock,'ddmmmyyyy_HHMM');
    s_name = ['emp_nom_rrl_greedy' s_date '.mat'];
    save(s_name,'-v7.3')    
    
    
end

%%  true cost
figure
all_epoch_cost = sum(epoch_costs(:,:,:), 2);
nom_epoch_cost = all_epoch_cost(1, 1, :);
rrl_epoch_cost = all_epoch_cost(2, 1, :);
greedy_epoch_cost = all_epoch_cost(3, 1, :);
boxplot([nom_epoch_cost(:) rrl_epoch_cost(:) greedy_epoch_cost(:)],'labels',{'nom','rrl', 'greedy'})
title('The true system cost')
grid on
%% worst case cost
figure
all_worst_cost = sum(worst_costs(:,:,:), 2);
nom_worst_cost = all_worst_cost(1, 1, :);
rrl_worst_cost = all_worst_cost(2, 1, :);
greedy_worst_cost = all_worst_cost(3, 1, :);
boxplot([nom_worst_cost(:) rrl_worst_cost(:) greedy_worst_cost(:)],'labels',{'nom','rrl', 'greedy'})
title('The worst-case cost')
grid on

%% uncertainty measure
figure
all_uncertainty_measure = sum(uncertainty_measure(:,:,:), 2);
nom_uncertainty_measure = all_uncertainty_measure(1, 1, :);
rrl_uncertainty_measure = all_uncertainty_measure(2, 1, :);
greedy_uncertainty_measure = all_uncertainty_measure(3, 1, :);
boxplot([nom_uncertainty_measure(:) rrl_uncertainty_measure(:) greedy_uncertainty_measure(:)],'labels',{'nom','rrl', 'greedy'})
title('The uncertainty measure')
grid on
