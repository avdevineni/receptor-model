% This script iterates receptor_model_v5 for varying bitter durations.

clear all

%% Model parameters

aM = 0.5; % closing rate (unbound or bound)
bM = 0; % opening rate (unbound)
bMs = 20; % opening rate (bound)
aH = 0.1; % inactivation rate (unbound)
aHs = 0.4; % inactivation rate (bound)
bH = 100; % de-inactivation rate (unbound)
bHs = 0; % de-inactivation rate (bound)

bitter_duration_list = [1 3 5 10];
num_expts = length(bitter_duration_list);

%% Run simulation

current_all = cell(1,num_expts);
m_all = cell(1,num_expts);
h_all = cell(1,num_expts);

for i = 1:num_expts
            
    bitter_duration = bitter_duration_list(i);
    
    [m,h,current] = receptor_model_main_func(bitter_duration,aM,bM,bMs,aH,aHs,bH,bHs);

    current_all{i} = current;
    m_all{i} = m;
    h_all{i} = h;            
end

%% Plot results

figure;
set(gcf, 'Position', [300,300,num_expts*350,250])

for i = 1:num_expts

        subplot(1,num_expts,i)
        plot(current_all{i},'LineWidth', 3, 'Color','k')
        maxval = 1.1;
        ylim([0 maxval])
        label1 = strcat('duration=',num2str(bitter_duration_list(i)));
        set(gca,'Xticklabel',[]) 
        title(label1,'FontSize',18)

end
    
% saveas(gcf,strcat('vary_duration','.fig')); 
% saveas(gcf,strcat('vary_duration','.png')); 
% 
% save('varyduration.mat', 'current_all','m_all','h_all',...
%     'bitter_duration_list','aM','bM','bMs','aH','aHs','bH','bHs')