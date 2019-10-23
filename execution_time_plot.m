% plot execution time

load('nom_rrl_greedy.mat');

%execution time of function robust controller
nom_data = execution_time(1,:,:);
%execution time of robust controller multistep run
nom_multistep = sum(execution_time(1,:,:), 2);
% execution time of multistep_rrl for time horizon 10
rrl_horizon = execution_time(2,1,:); 
 % execution time of multistep_rrl for time horizon 10 + execution time of
 % multipliers
rrl_horizon_multipliers = execution_time(2,1,:)+execution_time(4,1,:);
% full multistep calculation of rrl
rrl_horizon_multistep = sum(execution_time(2,:,:),2)+sum(execution_time(4,1,:),2);
% execution time of function optimal noise
optimal_noise =  execution_time(3,:,:);
% execution time of greedy (robust controller + optimal noise)
greedy_data =  execution_time(3,:,:) +  execution_time(5,:,:);
% full multistep calculation of greedy
greedy_multistep = sum(greedy_data, 2);

%%

fprintf('\nExecution time mean\n')
fprintf('\tnominal robust controller single = %.5e\n',mean(nom_data(:)))
fprintf('\tnominal robust controller multistep = %.5e\n',mean(nom_multistep(:)))
fprintf('\trrl single with time horizon 10   = %.5e\n',mean(rrl_horizon(:)))
fprintf('\trrl multiple run  = %.5e\n',mean(rrl_horizon_multistep(:)))
fprintf('\tgreedy single   = %.5e\n',mean(greedy_data(:)))
fprintf('\tgreedy multiple  = %.5e\n',mean(greedy_multistep(:)))

fprintf('\nExecution time std\n')
fprintf('\tnominal robust controller single = %.5e\n',std(nom_data(:)))
fprintf('\tnominal robust controller multistep = %.5e\n',std(nom_multistep(:)))
fprintf('\trrl single with time horizon 10   = %.5e\n',std(rrl_horizon(:)))
fprintf('\trrl multiple run  = %.5e\n',std(rrl_horizon_multistep(:)))
fprintf('\tgreedy single   = %.5e\n',std(greedy_data(:)))
fprintf('\tgreedy multiple  = %.5e\n',std(greedy_multistep(:)))

%% plot single and multiple run
figure
nom_data = nom_data(1,1,:);
greedy_data = greedy_data(1,1,:);
boxplot([nom_data(:) rrl_horizon(:) greedy_data(:)],'labels',{'nom','rrl', 'greedy'})
title('Single step')
grid on

figure
boxplot([nom_multistep(:) rrl_horizon_multistep(:) greedy_multistep(:)],'labels',{'nom','rrl', 'greedy'})
title('Multistep')
grid on

