%% Find Pupil, Blink, and Microsaccade Events For Blank/Jitter/ISI Trials - Sorted EyeLink Trials

% This code will look for certain EyeLink data types and try to recreate
% the Eyelink data dynamics found in PP and PnP trials. EyeLink timecourses
% are made based on these sorted trials. 

%Written by: Sharif I. Kronemer
%Date: 5/25/2022
%Last modified: 7/10/2022

clear

%% Prompts

%Select run location
%prompt_1 = 'Running code local or server [l, s]: ';
run_location = 's';%input(prompt_1,'s');

%Modality
prompt_3 = 'Modality [EEG, MRI]: ';
modality = input(prompt_3,'s');

% %EEG data is split between relevant/irrelevant conditions
% if isequal(modality, 'EEG')
    
    %Relevant condition
    rel_condition = input('Relevant or Irrelevant [r, i]: ','s');
    
% %MRI data is always rel
% else
%     
%     rel_condition = 'r';
% 
% end

%Location condition
prompt_2 = 'Center or Quadrant [c, q]: ';
condition = input(prompt_2,'s');

%Center relevant
if isequal(condition, 'c')

    irrelevant_location = 'Center Irrelevant';
    relevant_location = 'Center Relevant';
    con_name = 'center';
    
%Quadrant relevant
elseif isequal(condition, 'q')
    
    irrelevant_location = 'Quadrant Irrelevant';
    relevant_location = 'Quadrant Relevant';
    con_name = 'quadrant';
    
end

%Trial type [jitter, blank]
trial_type = 'jitter';

%%  Define directories and paths

%Save direcotry
if isequal(run_location,'l')

elseif isequal(run_location,'s')
    
    %Add paths
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/Machine Learning/Ensemble SVM Approach')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/Machine Learning')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/Behavioral Analysis')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/No Report Paradigm Machine Learning/Training Testing and Plotting Functions/MRI functions')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/No Report Paradigm Machine Learning/Training Testing and Plotting Functions/EEG functions')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/No Report Paradigm Machine Learning/Training Testing and Plotting Functions/Eyelink functions')

    %Save directory
    if strcmp(rel_condition, 'i')
        
        save_name = 'irrel';
            
    elseif strcmp(rel_condition, 'r')
        
        save_name = 'rel';
        
    end
    
    %Report group Eyelink data
    report_data_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Data/EyeLink';
    
    %Subject dir
    if isequal(modality, 'EEG')
        
        sub_dir = '/mnt/Data28/HNCT No Report Paradigm/Subject Analysis EEG';
    
    elseif isequal(modality, 'MRI')
        
        sub_dir = '/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI';
    
    end
    
    %EyeLink selected event save_dir
    save_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/EyeLink Analysis/EyeLink Event Times',modality,'Jitter Epochs');

end

%Make save directory
mkdir(save_dir)

%% STEP 1: Load EEG and Eyelink data - Combine over subjects

%Load subject EEG and EEG-eyelink data and add to a group variable. Output
%the group extracted eyelink, EEG, and true CP/CnP labels. 
if isequal(modality, 'EEG')
    
    disp('Extract EEG Eyelink Data')
   
    % Task Irrelevant
    if strcmp(rel_condition,'i')

        %Extract the irrelevant eyelink trials - EEG
        [noreport_data_pupil, noreport_data_tsfreshpupil, noreport_data_Ygaze, noreport_data_Xgaze, noreport_data_blink, noreport_data_microsac, group_subject,...
            group_run_number, group_trial_number, group_condition, group_relevant_location, group_relevant_trial, number_of_EEG_subjects, stimulus_trial_opacity, subject_trial_index] = ...
            combine_eyelink_EEG_irrelevant_epochs_jitter_function(sub_dir, irrelevant_location, run_location);
        
    % Task Relevant
    elseif strcmp(rel_condition,'r')
        
        %Extract the relevant eyelink trials - EEG
        [noreport_data_pupil, noreport_data_tsfreshpupil, noreport_data_Ygaze, noreport_data_Xgaze, noreport_data_blink, noreport_data_microsac, true_labels, group_subject,...
            group_run_number, group_trial_number, group_condition, group_relevant_location, group_relevant_trial, number_of_EEG_subjects, stimulus_trial_opacity, subject_trial_index] = ...
            combine_eyelink_EEG_relevant_epochs_jitter_function(sub_dir, relevant_location, run_location);
        
    end
    
elseif isequal(modality, 'MRI')
        
    disp('Extract MRI Eyelink Data')
        
    % Task Irrelevant
    if strcmp(rel_condition,'i')
    
        [noreport_data_pupil, noreport_data_tsfreshpupil, noreport_data_Ygaze, noreport_data_Xgaze, noreport_data_blink, noreport_data_microsac, group_subject,...
            group_run_number, group_trial_number, group_condition, group_relevant_location, group_BOLD_epoch_count, number_of_MRI_subjects, subject_trial_index, stimulus_trial_opacity] = ...
            combine_eyelink_MRI_irrelevant_epochs_function(sub_dir, irrelevant_location, run_location);
    
    % Task Relevant
    elseif strcmp(rel_condition,'r')

        [noreport_data_pupil, noreport_data_Ygaze, noreport_data_Xgaze, noreport_data_blink, noreport_data_microsac, group_subject, group_run_number, ... 
            group_trial_number, group_condition, group_relevant_location, group_relevant_trial, group_BOLD_epoch_count, number_of_MRI_subjects, subject_trial_index, stimulus_trial_opacity] = ...
            combine_eyelink_MRI_rel_irrel_epochs_jitter_function(sub_dir, relevant_location, run_location);
        
    end
    
end   

%% Pupil Epochs

%Initialize variables
big_pupil_trials = [];
small_pupil_trials = [];
pupil_trial_condition_idx = [];

%Distance from center of trial
if isequal(trial_type, 'blank')
    
    %Search window
    onset_window = 7001;
    offset_window = 8000;
    
    % Baseline 
    baseline_window = 5001:6000;

elseif isequal(trial_type, 'jitter')

    %Search window
    onset_window = 10001; 
    offset_window = 11000; 
    
    % Baseline
    baseline_window = 8001:9000; 
    
end

%Number of samples
sample_num = 25;

%Loop over blink trials
for trial = 1:size(noreport_data_pupil,1)
    
    %Calculate the SD/mean values
    trial_SD = std(noreport_data_pupil(trial,:));
    trial_mean = nanmean(noreport_data_pupil(trial,:));
    
    if sum(noreport_data_pupil(trial,onset_window:offset_window) > trial_mean + trial_SD/3) > sample_num && ...
        sum(noreport_data_pupil(trial,onset_window-500:onset_window) < trial_mean + trial_SD/3) > sample_num
            
         %Add idx
        pupil_trial_condition_idx(trial,1) = 1;
    
        
    else
        
        %Add idx
        pupil_trial_condition_idx(trial,1) = 0;
    
    end    

end

% %Randomly mix some percentage of trials
% num_trials = 0;
% rand_idx = find(pupil_trial_condition_idx == 1);
% rand_idx = sort(rand_idx(randperm(length(rand_idx),num_trials)));
% 
% %Flip index
% pupil_trial_condition_idx(rand_idx) = 0;

% Extract trials
pupil_pupil_sorted_trials = noreport_data_pupil(pupil_trial_condition_idx == 1,:);
no_pupil_pupil_sorted_trials = noreport_data_pupil(pupil_trial_condition_idx == 0,:);

blink_pupil_sorted_trials = noreport_data_blink(pupil_trial_condition_idx == 1,:);
no_blink_pupil_sorted_trials = noreport_data_blink(pupil_trial_condition_idx == 0,:);

microsac_pupil_sorted_trials = noreport_data_microsac(pupil_trial_condition_idx == 1,:);
no_microsac_pupil_sorted_trials = noreport_data_microsac(pupil_trial_condition_idx == 0,:);

% big_pupil_trials = noreport_data_pupil(pupil_trial_condition_idx == 1,:);
% small_pupil_trials = noreport_data_pupil(pupil_trial_condition_idx == 0,:);
% 
% blink_trials = noreport_data_blink(pupil_trial_condition_idx == 1,:);
% no_blink_trials = noreport_data_blink(pupil_trial_condition_idx == 0,:);
% 
% microsaccade_trials = noreport_data_microsac(pupil_trial_condition_idx == 1,:);
% no_microsaccade_trials = noreport_data_microsac(pupil_trial_condition_idx == 0,:);

%{
%Create a baseline matrix across subjects
big_baseline_matrix = nanmean(big_pupil_trials(:,baseline_window),2);
small_baseline_matrix = nanmean(small_pupil_trials(:,baseline_window),2);

%Baseline
big_pupil_baselined = big_pupil_trials - shiftdim(repmat(big_baseline_matrix,1,12000),1)';
small_pupil_baselined = small_pupil_trials - shiftdim(repmat(small_baseline_matrix,1,12000),1)';

% Smooth data
bin_size = 500;

smoothed_microsaccade_trials = movmean(microsaccade_trials',bin_size)';
smoothed_no_microsaccade_trials = movmean(no_microsaccade_trials',bin_size)';

% Plot mean timecourses - PUPIL
pupil_fig = figure
hold on

title(['Pupil Sorted - Pupil Timecourse ',modality, ' ',con_name,' ',save_name])
xlabel('Time (ms)')
ylabel('Baselined Diameter (mm)')

ylim([-0.15 0.15])
%xlim([-1000 5000])
xlim([-1000 6000])

yticks([-0.15:0.05:0.15])

timevector = [-1000:5999];%[-1000:4999];
time_period = [5001:12000];%[5001:11000];

plot([0 0],[-0.15 0.15],'k')

big = plot(timevector, nanmean(big_pupil_baselined(:,time_period),1),'b')
small = plot(timevector, nanmean(small_pupil_baselined(:,time_period),1),'r')

legend([big,small],{['Big N = ',num2str(size(big_pupil_trials,1))], ['Small N = ',num2str(size(small_pupil_trials,1))]})

%Save figure
cd(save_dir)
savefig(pupil_fig,['Pupil_sorted_pupil_timecourse_',con_name,'_',save_name,'_selected_trials.fig'])

% Plot mean timecourses - BLINK
blink_fig = figure
hold on

title(['Pupil Sorted - Blink Timecourse ',modality, ' ',con_name,' ',save_name])
xlabel('Time (ms)')
ylabel('Blink occurrence')

ylim([0 0.3])
%xlim([-1000 5000])
xlim([-1000 6000])

yticks([0.02:0.04:0.3])

timevector = [-1000:5999];%[-1000:4999];
time_period = [5001:12000];%[5001:11000];

plot([0 0],[0,0.3],'k')

blink = plot(timevector, nanmean(blink_trials(:,time_period),1), 'b')
no_blink = plot(timevector, nanmean(no_blink_trials(:,time_period),1), 'r')

legend([blink,no_blink],{['Blink N = ',num2str(size(blink_trials,1))], ['No Blink N = ',num2str(size(no_blink_trials,1))]})

%Save figure
cd(save_dir)
savefig(blink_fig,['Pupil_sorted_blink_timecourse_',con_name,'_',save_name,'_selected_trials.fig'])

% Plot mean timecourses - MICROSACCADE
micro_fig = figure
hold on

title(['Pupil Sorted - Microsaccade Timecourse ',modality, ' ',con_name,' ',save_name])
xlabel('Time (ms)')
ylabel('Microsaccade (%)')

yticks([0.002:0.002:0.018])
ylim([0.002 0.018])

timevector = [-1000:5999];%[-1000:4999];
time_period = [5001:12000];%[5001:11000];

plot([0,0],[0,0.018],'k')

micro = plot(timevector, nanmean(smoothed_microsaccade_trials(:,time_period),1),'r')
no_micro = plot(timevector, nanmean(smoothed_no_microsaccade_trials(:,time_period),1),'b')

legend([micro,no_micro],{['Microsaccade N = ',num2str(size(microsaccade_trials,1))], ['No Microsaccade N = ',num2str(size(no_microsaccade_trials,1))]})

%Save figure
cd(save_dir)
savefig(micro_fig,['Pupil_sorted_microsaccade_timecourse_',con_name,'_',save_name,'_selected_trials.fig'])

close all
%}

%% Blink Epochs

%Distance from center of trial
if isequal(trial_type, 'blank')
    
    %Search window
    onset_window = 7001;
    offset_window = 8000;
    
    % Baseline 
    baseline_window = 5001:6000;

elseif isequal(trial_type, 'jitter')

    %Search window
    onset_window = 10001; 
    offset_window = 11000; 
    
    % Baseline
    baseline_window = 8001:9000; 
    
end

%Number of samples
sample_num = 25;

%Initialize variables
blink_trial_condition_idx = [];

%Loop over blink trials
for trial = 1:size(noreport_data_blink,1)
    
    % Blink present
    if sum(noreport_data_blink(trial,onset_window:offset_window)) > sample_num
               
        %Add idx
        blink_trial_condition_idx(trial,1) = 1;
    
    % Blink absent
    else
     
        % Add idx
        blink_trial_condition_idx(trial,1) = 0;
        
    end    

end

% %Randomly mix some percentage of trials
% rand_idx = find(blink_trial_condition_idx == 1);
% num_trials = round(length(rand_idx)*.25);
% rand_idx = sort(rand_idx(randperm(length(rand_idx),num_trials)));
% 
% %Flip index
% blink_trial_condition_idx(rand_idx) = 0;

% Extract trials
pupil_blink_sorted_trials = noreport_data_pupil(blink_trial_condition_idx == 1,:);
no_pupil_blink_sorted_trials = noreport_data_pupil(blink_trial_condition_idx == 0,:);

blink_blink_sorted_trials = noreport_data_blink(blink_trial_condition_idx == 1,:);
no_blink_blink_sorted_trials = noreport_data_blink(blink_trial_condition_idx == 0,:);

microsac_blink_sorted_trials = noreport_data_microsac(blink_trial_condition_idx == 1,:);
no_microsac_blink_sorted_trials = noreport_data_microsac(blink_trial_condition_idx == 0,:);

% big_pupil_trials = noreport_data_pupil(blink_trial_condition_idx == 1,:);
% small_pupil_trials = noreport_data_pupil(blink_trial_condition_idx == 0,:);
% 
% blink_trials = noreport_data_blink(blink_trial_condition_idx == 1,:);
% no_blink_trials = noreport_data_blink(blink_trial_condition_idx == 0,:);
% 
% microsaccade_trials = noreport_data_microsac(blink_trial_condition_idx == 1,:);
% no_microsaccade_trials = noreport_data_microsac(blink_trial_condition_idx == 0,:);

%{
%Create a baseline matrix across subjects
big_baseline_matrix = nanmean(big_pupil_trials(:,baseline_window),2);
small_baseline_matrix = nanmean(small_pupil_trials(:,baseline_window),2);

%Baseline
big_pupil_baselined = big_pupil_trials - shiftdim(repmat(big_baseline_matrix,1,12000),1)';
small_pupil_baselined = small_pupil_trials - shiftdim(repmat(small_baseline_matrix,1,12000),1)';

% Smooth data
bin_size = 500;

smoothed_microsaccade_trials = movmean(microsaccade_trials',bin_size)';
smoothed_no_microsaccade_trials = movmean(no_microsaccade_trials',bin_size)';

% Plot mean timecourses - PUPIL
pupil_fig = figure
hold on

title(['Blink Sorted - Pupil Timecourse ',modality, ' ',con_name,' ',save_name])
xlabel('Time (ms)')
ylabel('Baselined Diameter (mm)')

ylim([-0.15 0.15])
%xlim([-1000 5000])
xlim([-1000 6000])

yticks([-0.15:0.05:0.15])

timevector = [-1000:5999];%[-1000:4999];
time_period = [5001:12000];%[5001:11000];

plot([0 0],[-0.15 0.15],'k')

big = plot(timevector, nanmean(big_pupil_baselined(:,time_period),1),'b')
small = plot(timevector, nanmean(small_pupil_baselined(:,time_period),1),'r')

legend([big,small],{['Big N = ',num2str(size(big_pupil_trials,1))], ['Small N = ',num2str(size(small_pupil_trials,1))]})

%Save figure
cd(save_dir)
savefig(pupil_fig,['Blink_sorted_pupil_timecourse_',con_name,'_',save_name,'_selected_trials.fig'])

% Plot mean timecourses - BLINK
blink_fig = figure
hold on

title(['Blink Sorted - Blink Timecourse ',modality, ' ',con_name,' ',save_name])
xlabel('Time (ms)')
ylabel('Blink occurrence')

ylim([0 0.3])
%xlim([-1000 5000])
xlim([-1000 6000])

yticks([0.02:0.04:0.3])

timevector = [-1000:5999];%[-1000:4999];
time_period = [5001:12000];%[5001:11000];

plot([0 0],[0,0.3],'k')

blink = plot(timevector, nanmean(blink_trials(:,time_period),1), 'b')
no_blink = plot(timevector, nanmean(no_blink_trials(:,time_period),1), 'r')

legend([blink,no_blink],{['Blink N = ',num2str(size(blink_trials,1))], ['No Blink N = ',num2str(size(no_blink_trials,1))]})

%Save figure
cd(save_dir)
savefig(blink_fig,['Blink_sorted_blink_timecourse_',con_name,'_',save_name,'_selected_trials.fig'])

% Plot mean timecourses - MICROSACCADE
micro_fig = figure
hold on

title(['Blink Sorted - Microsaccade Timecourse ',modality, ' ',con_name,' ',save_name])
xlabel('Time (ms)')
ylabel('Microsaccade (%)')

yticks([0.002:0.002:0.018])
ylim([0.002 0.018])

timevector = [-1000:5999];%[-1000:4999];
time_period = [5001:12000];%[5001:11000];

plot([0,0],[0,0.018],'k')

micro = plot(timevector, nanmean(smoothed_microsaccade_trials(:,time_period),1),'r')
no_micro = plot(timevector, nanmean(smoothed_no_microsaccade_trials(:,time_period),1),'b')

legend([micro,no_micro],{['Microsaccade N = ',num2str(size(microsaccade_trials,1))], ['No Microsaccade N = ',num2str(size(no_microsaccade_trials,1))]})

%Save figure
cd(save_dir)
savefig(micro_fig,['Blink_sorted_microsaccade_timecourse_',con_name,'_',save_name,'_selected_trials.fig'])
%}

%% Microsaccade Epochs

%Distance from center of trial
if isequal(trial_type, 'blank')
    
    %Search window
    onset_window = 6251;
    offset_window = 6875;
    
    % Baseline 
    baseline_window = 5001:6000;

elseif isequal(trial_type, 'jitter')

    %Search window
    onset_window = 9251; 
    offset_window = 9875; 
    
    % Baseline
    baseline_window = 8001:9000; 
    
end

%Initialize variables
%microsaccade_trials = [];
%no_microsaccade_trials = [];
microsaccade_trial_condition_idx = [];

%Number of samples
sample_num = 5;

%Loop over blink trials
for trial = 1:size(noreport_data_microsac,1)
    
    if sum(noreport_data_microsac(trial,onset_window:offset_window)) < sample_num 
                
        %Add idx
        microsaccade_trial_condition_idx(trial,1) = 1;
    
    else
                
        %Add idx
        microsaccade_trial_condition_idx(trial,1) = 0;
    
        
    end    

end

% %Randomly mix some percentage of trials
% post_rand_idx = find(microsaccade_trial_condition_idx == 1);
% neg_rand_idx = find(microsaccade_trial_condition_idx == 0);
% 
% num_trials = round(length(post_rand_idx)*0.25);
% post_rand_idx = sort(post_rand_idx(randperm(length(post_rand_idx),num_trials)));
% neg_rand_idx = sort(neg_rand_idx(randperm(length(neg_rand_idx),num_trials)));
% 
% %Flip index
% microsaccade_trial_condition_idx(post_rand_idx) = 0;
% microsaccade_trial_condition_idx(neg_rand_idx) = 1;

% Extract trials
pupil_microsac_sorted_trials = noreport_data_pupil(microsaccade_trial_condition_idx == 1,:);
no_pupil_microsac_sorted_trials = noreport_data_pupil(microsaccade_trial_condition_idx == 0,:);

blink_microsac_sorted_trials = noreport_data_blink(microsaccade_trial_condition_idx == 1,:);
no_blink_microsac_sorted_trials = noreport_data_blink(microsaccade_trial_condition_idx == 0,:);

microsac_microsac_sorted_trials = noreport_data_microsac(microsaccade_trial_condition_idx == 1,:);
no_microsac_microsac_sorted_trials = noreport_data_microsac(microsaccade_trial_condition_idx == 0,:);

% big_pupil_trials = noreport_data_pupil(microsaccade_trial_condition_idx == 1,:);
% small_pupil_trials = noreport_data_pupil(microsaccade_trial_condition_idx == 0,:);
% 
% blink_trials = noreport_data_blink(microsaccade_trial_condition_idx == 1,:);
% no_blink_trials = noreport_data_blink(microsaccade_trial_condition_idx == 0,:);
% 
% microsaccade_trials = noreport_data_microsac(microsaccade_trial_condition_idx == 0,:);
% no_microsaccade_trials = noreport_data_microsac(microsaccade_trial_condition_idx == 1,:);

%{
%Create a baseline matrix across subjects
big_baseline_matrix = nanmean(big_pupil_trials(:,baseline_window),2);
small_baseline_matrix = nanmean(small_pupil_trials(:,baseline_window),2);

%Baseline
big_pupil_baselined = big_pupil_trials - shiftdim(repmat(big_baseline_matrix,1,12000),1)';
small_pupil_baselined = small_pupil_trials - shiftdim(repmat(small_baseline_matrix,1,12000),1)';

% Smooth data
bin_size = 500;

smoothed_microsaccade_trials = movmean(microsaccade_trials',bin_size)';
smoothed_no_microsaccade_trials = movmean(no_microsaccade_trials',bin_size)';

% Plot mean timecourses - PUPIL
pupil_fig = figure
hold on

title(['Microsaccade Sorted - Pupil Timecourse ',modality, ' ',con_name,' ',save_name])
xlabel('Time (ms)')
ylabel('Baselined Diameter (mm)')

ylim([-0.15 0.15])
%xlim([-1000 5000])
xlim([-1000 6000])

yticks([-0.15:0.05:0.15])

timevector = [-1000:5999];%[-1000:4999];
time_period = [5001:12000];%[5001:11000];

plot([0 0],[-0.15 0.15],'k')

big = plot(timevector, nanmean(big_pupil_baselined(:,time_period),1),'b')
small = plot(timevector, nanmean(small_pupil_baselined(:,time_period),1),'r')

legend([big,small],{['Big N = ',num2str(size(big_pupil_trials,1))], ['Small N = ',num2str(size(small_pupil_trials,1))]})

%Save figure
cd(save_dir)
savefig(pupil_fig,['Microsaccade_sorted_pupil_timecourse_',con_name,'_',save_name,'_selected_trials.fig'])

% Plot mean timecourses - BLINK
blink_fig = figure
hold on

title(['Microsaccade Sorted - Blink Timecourse ',modality, ' ',con_name,' ',save_name])
xlabel('Time (ms)')
ylabel('Blink occurrence')

ylim([0 0.3])
%xlim([-1000 5000])
xlim([-1000 6000])

yticks([0.02:0.04:0.3])

timevector = [-1000:5999];%[-1000:4999];
time_period = [5001:12000];%[5001:11000];

plot([0 0],[0,0.3],'k')

blink = plot(timevector, nanmean(blink_trials(:,time_period),1), 'b')
no_blink = plot(timevector, nanmean(no_blink_trials(:,time_period),1), 'r')

legend([blink,no_blink],{['Blink N = ',num2str(size(blink_trials,1))], ['No Blink N = ',num2str(size(no_blink_trials,1))]})

%Save figure
cd(save_dir)
savefig(blink_fig,['Microsaccade_sorted_blink_timecourse_',con_name,'_',save_name,'_selected_trials.fig'])

% Plot mean timecourses - MICROSACCADE
micro_fig = figure
hold on

title(['Microsaccade Sorted - Microsaccade Timecourse ',modality, ' ',con_name,' ',save_name])
xlabel('Time (ms)')
ylabel('Microsaccade (%)')

yticks([0.002:0.002:0.018])
ylim([0.002 0.018])

timevector = [-1000:5999];%[-1000:4999];
time_period = [5001:12000];%[5001:11000];

plot([0,0],[0,0.018],'k')

micro = plot(timevector, nanmean(smoothed_microsaccade_trials(:,time_period),1),'r')
no_micro = plot(timevector, nanmean(smoothed_no_microsaccade_trials(:,time_period),1),'b')

legend([micro,no_micro],{['Microsaccade N = ',num2str(size(microsaccade_trials,1))], ['No Microsaccade N = ',num2str(size(no_microsaccade_trials,1))]})

%Save figure
cd(save_dir)
savefig(micro_fig,['Microsaccade_sorted_microsaccade_timecourse_',con_name,'_',save_name,'_selected_trials.fig'])
%}

%% Save Data

cd(save_dir)
    
if isequal(trial_type, 'blank')
    
    save(['eyelink_',con_name,'_',save_name,'_selected_blank_trials.mat'],  'pupil_trial_condition*', 'blink_trial_condition*', ...
        'microsaccade_trial_condition*', 'subject_trial_index', 'stimulus_trial_opacity', 'noreport*')
    
elseif isequal(trial_type, 'jitter')
    
    
    % Task Relevant
    if strcmp(rel_condition,'r')
        
        save(['eyelink_',con_name,'_',save_name,'_selected_jitter_trials.mat'], 'pupil_trial_condition*', 'blink_trial_condition*', 'group_condition', 'group_relevant_trial', ...
            'microsaccade_trial_condition*', 'subject_trial_index', 'stimulus_trial_opacity', 'noreport*', 'pupil*', 'blink*', 'microsac*', 'no_pupil*', 'no_blink*', 'no_microsac*')

    % Task Irrelevant
    elseif strcmp(rel_condition,'i')
        
        save(['eyelink_',con_name,'_',save_name,'_selected_jitter_trials.mat'], 'pupil_trial_condition*', 'blink_trial_condition*', 'group_condition', ...
            'microsaccade_trial_condition*', 'subject_trial_index', 'stimulus_trial_opacity', 'noreport*', 'pupil*', 'blink*', 'microsac*', 'no_pupil*', 'no_blink*', 'no_microsac*')
            
    end

end
