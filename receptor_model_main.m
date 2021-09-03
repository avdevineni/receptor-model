% This is a single receptor model where open receptors gradually become
% inactivated, and inactivation is relieved by bitter removal. 

% This version is a deterministic model. Current is given by the product of
% variable m (reflects receptor opening) and variable h (reflects
% receptor's inactivation state). This models the behavior of a population
% of receptors.

% Receptors can be open and active (m=1,h=1), open and inactive (m=1,h=0),
% closed and active (m=0,h=1), or closed and inactive (m=0,h=0).

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

%% Plot results

figure
set(gcf, 'Position', [300,400,400,700])

subplot(4,1,1)
plot(t,s,'k','linewidth',2)
ylim([0 1.1])
title('Ligand binding')

subplot(4,1,2)
plot(t,m,'k','linewidth',2)
ylim([0 1.1])
title('Open receptors')

subplot(4,1,3)
plot(t,1-h,'k','linewidth',2)
ylim([0 1.1])
title('Inactivated receptors')

subplot(4,1,4)
plot(t,current,'k','linewidth',2)
ylim([0 1.1])
title('Receptor current')
xlabel('time (s)')

% saveas(gcf,'main_results.fig')
% saveas(gcf,'main_results.png')
% 
% save('results.mat', 'current','m','h',...
%     'bitter_duration','aM','bM','bMs','aH','aHs','bH','bHs')