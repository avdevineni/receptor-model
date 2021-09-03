function [m,h,current] = receptor_model_main_func(bitter_duration,aM,bM,bMs,aH,aHs,bH,bHs)

% Function to iterate receptor_model_v5

T = bitter_duration + 15; % total time in sec
dt = 0.001; % size of timestep
t = 0:dt:T;

m = zeros(length(t),1); % all receptors start closed
h = ones(length(t),1); % all receptors start in non-inactivated state

s = double(((t>5) & (t<= (5+bitter_duration)))); % bitter presented at t=5 sec

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

end

