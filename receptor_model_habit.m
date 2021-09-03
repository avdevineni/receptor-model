% This is a single receptor model where open receptors gradually become
% inactivated, and inactivation is relieved by bitter removal. 

% This version is a deterministic model. Current is given by the product of
% variable m (reflects receptor opening) and variable h (reflects
% receptor's inactivation state). This models the behavior of a population
% of receptors.

% Receptors can be open and active (m=1,h=1), open and inactive (m=1,h=0),
% closed and active (m=0,h=1), or closed and inactive (m=0,h=0).

% This script simulates repeated bitter stimulation to model habituation by
% modulating current based on the number of previous stimuli and the
% receptor's bound state.

%% Model parameters

aM = 0.5; % closing rate (unbound or bound)
bM = 0; % opening rate (unbound)
bMs = 20; % opening rate (bound)
aH = 0.1; % inactivation rate (unbound)
aHs = 0.4; % inactivation rate (bound)
bH = 100; % de-inactivation rate (unbound)
bHs = 0; % de-inactivation rate (bound)

bitter_duration = 5; % in sec
T = bitter_duration + 15; % total time in sec
dt = 0.001; % size of timestep
t = 0:dt:T;

m = zeros(length(t),1); % all receptors start closed
h = ones(length(t),1); % all receptors start in non-inactivated state

s = double(((t>5) & (t<= (5+bitter_duration)))); % bitter presented at t=5 sec

% habituation parameters
num_stim = 10;
tau_H_bound = 1.8; % tau for response decay if bound
tau_H_unbound = 15; % tau for response decay if unbound
H_bound = zeros(1,num_stim);
H_unbound = zeros(1,num_stim);
for i = 1:num_stim
    H_bound(i) = exp(-(i-1)./tau_H_bound);
    H_unbound(i) = exp(-(i-1)./tau_H_unbound);
end

% H determines how current is modulated by habituation
H = zeros(num_stim, length(t));
for n = 1:num_stim
    for i = 1:length(t)
        if s(i) == 0
            H(n,i) = H_unbound(n);
        else
            H(n,i) = H_bound(n);
        end
    end
end


%% Run simulation

for i=2:length(t)
    
    % choose which values to use depending on if bound or unbound
    if s(i) == 1
        bM_curr = bMs;
        aH_curr = aHs;
        bH_curr = bHs;
    else
        bM_curr = bM;
        aH_curr = aH;
        bH_curr = bH;
    end
        
    % m(t) is calculated as m(t-1), minus open receptors that closed, plus
    % closed receptors that opened
    m(i) = m(i-1) - aM*dt*m(i-1) + bM_curr*dt*(1-m(i-1)); 
    
    % h(t) is calculated as h(t-1), minus active receptors that
    % inactivated, plus inactive receptors that de-inactivated
    h(i) = h(i-1) - aH_curr*dt*h(i-1) + bH_curr*dt*(1-h(i-1)); 
    
end

current = m.*h;
current_alltrials = repmat(current',num_stim,1);
current_habit = current_alltrials .* H;

%% Plot results

figure;
set(gcf, 'Position', [100,400,1200,300])

subplot(1,3,1)
plot(t,m,'LineWidth', 2, 'Color','r')
hold on
plot(t,1-h,'LineWidth', 2, 'Color','b')
ylim([0 1.1])
xlim([0 T])
title('Receptor states','FontSize',18)
legend('open receptors', 'inactive receptors')
xlabel('time (s)')

subplot(1,3,2)
plot(t,current,'LineWidth', 2, 'Color','k')
ylim([0 1.1])
title('Receptor current','FontSize',18)
xlabel('time (s)')

subplot(1,3,3)
plot(t,current_habit','LineWidth', 2)
ylim([0 1.1])
title('Response to repeated stim','FontSize',18)
legend(num2str(linspace(1,num_stim,num_stim)')) 
xlabel('time (s)')

% saveas(gcf,'habit.fig')
% saveas(gcf,'habit.png')
% 
% save('habit.mat', 'current','current_habit','m','h',...
%     'bitter_duration','aM','bM','bMs','aH','aHs','bH','bHs')