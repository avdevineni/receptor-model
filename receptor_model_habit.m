% This is a single receptor model where open receptors gradually become
% inactivated, and inactivation is terminated by bitter removal. 
% The prolonged OFF response is due to receptors closing slowly.

% This version is modified to add a process that could explain differential
% habituation of ON vs OFF response. Assume there is a process that
% suppresses the output of the receptor but does so differentially depending
% on whether bitter is bound. This is implemented globally, not at the
% individual receptor level.

clear all

%% Set parameters

% Parameters 
rng(0); % sets random seed
num_receptors = 10000; % number of receptors
baseline_open = 0; % initial fraction of receptors open (at baseline)

dt = 0.01; % size of timestep

% Probabilities of channels opening and closing at each timestep
prob_open_baseline = 0; % probability of closed channel opening (unbound)
prob_close_baseline = 0.002; % probability of open channel closing (unbound)
prob_open_bitter = 0.1; % probability of closed channel opening (when bound)
prob_close_bitter = 0.002; % probability of open channel closing (when bound)
prob_inact_baseline = 0.001; % probability of open channel inactivating (unbound)
prob_inact_bitter = 0.006; % probability of open channel inactivating (bound)
prob_deinact = 0.0001; % probability of open channel de-inactivating (same for bound or unbound)
prob_deinact_removal = 0.9; % probability of open channel de-inactivating (upon bitter removal) 

% arrays with probabilities for different epochs
open_probs = [prob_open_baseline prob_open_bitter];
close_probs = [prob_close_baseline prob_close_bitter];
inact_probs = [prob_inact_baseline prob_inact_bitter];
deinact_probs = [prob_deinact prob_deinact_removal];

% set up timing of simulation: 20 sec with 5 sec bitter at t=5 sec
bitter_starttime = 5; % in sec
bitter_duration = 5; % in sec
total_time = bitter_starttime + bitter_duration + 10; % in sec

num_timesteps = total_time/dt; % number of timesteps to iterate through

% set up bitter presentation 
bitter_startframe = bitter_starttime./dt + 1; % start on next frame
bitter_endframe = (bitter_starttime + bitter_duration)/dt;
bitter_presentation = zeros(num_timesteps,1); % determines when bitter is presented
bitter_presentation(bitter_startframe:bitter_endframe) = ones(1,bitter_duration./dt);

% decaying exponential to mimic bitter removal
expcoeff = 0.02; % tau in sec (determines rate of decay)
xvals = linspace(0, num_timesteps - bitter_endframe, num_timesteps - bitter_endframe + 1);
yvals = exp(-xvals ./ expcoeff .* dt);
yvals2 = yvals./sum(yvals); % want the vector to sum to 1

% spread out bitter removal from receptors according to this curve
bitter_removal = round(yvals2.*num_receptors);
bitter_removal(1) = bitter_removal(1) + (num_receptors - sum(bitter_removal)); % add extras to first frame
bitter_removal = nonzeros(bitter_removal); % remove zeros at end 

% list of which frame (starting from bitter endframe) each receptor has
% bitter removed
removal_times = [];
i = 1;
ind = 1;
while i < (num_receptors+1)
    removal_times = [removal_times ind*ones(1,bitter_removal(ind))];
    i = i + bitter_removal(ind);
    ind = ind+1;
end

removal_times = removal_times + bitter_endframe; % convert to frames

%% Run simulation

% Arrays to store state of each receptor at each time
receptor_openstates = zeros(num_receptors, num_timesteps); % open state
receptor_actstates = zeros(num_receptors, num_timesteps); % activation state

% Initialize receptor states (not really needed if all receptors start closed)
% 0 = closed, 1 = open
receptors = [zeros(1,(1-baseline_open).*num_receptors) ones(1,(baseline_open.*num_receptors))];  % Fill the vector with 0 and 1
receptors = receptors(randperm(num_receptors)); % Permute the zeros and ones

receptor_openstates(:,1) = receptors;

% Start simulation
for i = 2:num_timesteps
    
    % generate random numbers for each receptor
    randnums_open = rand(num_receptors, 1);
    randnums_inact = rand(num_receptors, 1);
    
    % figure out which epoch we're in to use appropriate probabilities
    if bitter_presentation(i) == 1 % during bitter
        prob_ind = 2;
    else % pre- or post-bitter
        prob_ind = 1;
    end
    
    for j = 1:num_receptors % iterate through all receptors
            
        % for closed receptors
        if receptor_openstates(j,i-1) == 0 % if closed
            if randnums_open(j) < open_probs(prob_ind)
                receptor_openstates(j,i) = 1; % open if needed
            end
        end

        % for open receptors
        if receptor_openstates(j,i-1) == 1 % if open

            if randnums_open(j) < close_probs(prob_ind)
                receptor_openstates(j,i) = 0; % close if needed

            else 
                receptor_openstates(j,i) = 1; % otherwise leave open

                % if leaving open: assess activation state and change if needed
                if receptor_actstates(j,i-1) == 0 % if not inactivated (=0)
                    if randnums_inact(j) < inact_probs(prob_ind)
                        receptor_actstates(j,i) = 1; % inactivate if needed (=1)
                    end

                else % if inactivated (=1)
                    
                    % if bitter was removed on this frame
                    if i == removal_times(j)
                    
                        % determine if receptor inactivation state changes
                        if randnums_inact(j) < deinact_probs(2) 
                            receptor_actstates(j,i) = 0; % de-inactivate (=0)
                        else
                            receptor_actstates(j,i) = 1; % keep inactivated (=1)
                        end
                        
                    else % if bitter not just removed
                        
                        if randnums_inact(j) < deinact_probs(1) 
                            receptor_actstates(j,i) = 0; % de-inactivate (=0)
                        else
                            receptor_actstates(j,i) = 1; % keep inactivated (=1)
                        end

                    end
                    
                end

            end
        end
            
    end
    
end

%% Calculate receptor activity before habituation

open_receptors = mean(receptor_openstates);
inact_receptors = mean(receptor_actstates);

act_states_new = (receptor_actstates - 1).*(-1); % flip 0s and 1s
receptor_activity = mean(receptor_openstates.*act_states_new);


%% Implement habituation
% Suppresses receptor output in a graded fashion, with different amounts of
% suppression depending on if receptor is bound/unbound

% Scaling factors determining how much to decrease output based on the
% number of previous bitter stimuli. 
% Model this with exponential decay.

% tau values set the rate of decay - different for ON and OFF
tau_habit_ON = 1.8; % approximated based on num of stim needed to decrease to 37% of max
tau_habit_OFF = 15; % approximated based on num of stim needed to decrease to 37% of max

xvals = linspace(0,1999,2000); 
ONhabit = exp(-xvals./tau_habit_ON * dt);
OFFhabit = exp(-xvals./tau_habit_OFF * dt);

num_stim = 10;

bound_scaling_factor = zeros(1,num_stim);
unbound_scaling_factor = zeros(1,num_stim);
for i = 1:num_stim
    xval = (i-1)*100+1;
    bound_scaling_factor(i) = ONhabit(xval);
    unbound_scaling_factor(i) = OFFhabit(xval);
end

scaling = zeros(num_stim, num_timesteps);
for n = 1:num_stim
    for i = 1:num_timesteps
        if bitter_presentation(i) == 0
            scaling(n,i) = unbound_scaling_factor(n);
        else
            scaling(n,i) = bound_scaling_factor(n);
        end
    end
end

traces = repmat(receptor_activity,num_stim,1);
final_output = traces .* scaling;

%% Plot results

figure;
set(gcf, 'Position', [100,400,1200,300])

subplot(1,3,1)
plot(open_receptors,'LineWidth', 2, 'Color','r')
hold on
plot(inact_receptors,'LineWidth', 2, 'Color','b')
ylim([0 1.1])
title('Receptor states','FontSize',18)
legend('open receptors', 'inactivated receptors')
set(gca,'Xticklabel',[]) 

subplot(1,3,2)
plot(receptor_activity,'LineWidth', 2, 'Color','k')
ylim([0 1.1])
title('Receptor activity','FontSize',18)
set(gca,'Xticklabel',[]) 

subplot(1,3,3)
plot(final_output','LineWidth', 2)
ylim([0 1.1])
title('Response to repeated stim','FontSize',18)
legend(num2str(linspace(1,num_stim,num_stim)'))
set(gca,'Xticklabel',[]) 

%saveas(gcf,'final_results_habit.fig')
%saveas(gcf,'final_results_habit.png')

save('habit_v3_2.mat', 'open_receptors','inact_receptors','receptor_activity','final_output')
