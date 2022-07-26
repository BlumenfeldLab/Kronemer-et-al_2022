%% Create Group Data of Eyelink Data - Report and No Report Paradigm PP and PnP Epochs

%This code will gather all the PP and PNP trials across subjects, average
%within subjects, and then create a group matrix of data type x time x
%subjects. The code also rejects subjects by poor behavioral and rejects
%trials with a blink event during the stimulus time. 

%This code is flexible for PP and PnP trials from task relevant and
%irrelevant conditions.

%Written by: Sharif I. Kronemer
%Date: 5/30/2021
%Last Modified: 5/15/2022

clear

%% Run Location

%Select run location
prompt_1 = 'Running code local or server [l, s]: ';
run_location = input(prompt_1,'s');

%Run Task Relevant or Irrrelevant Stimuli
prompt_2 = 'Relevant or irrelevant stimuli [rel, irrel, combine]: ';
relevant_condition = input(prompt_2, 's');

%Stimulus location set
prompt_3 = 'Stimulus location set [c, q, combine]: ';
location_set = input(prompt_3, 's');

%Stimulus Opacity
prompt_5 = 'Stimulus opacity (threshold, blank): ';
opacity = input(prompt_5, 's');

%Setup condition identifiers

%Relevant stimuli
if isequal(relevant_condition, 'rel')
    
    %Condition name
    condition_name = 'relevant';
    
    %Center Relevant
    if isequal(location_set, 'c')

        location_cell = {'Center Relevant'};
        rel_save_name = 'cent';

    %Quadrant Relevant
    elseif isequal(location_set, 'q')

        location_cell = {'Quadrant Relevant'};
        rel_save_name = 'quad';

    %Combine Center and Quadrant Relevant   
    elseif isequal(location_set, 'combine')

        location_cell = {'Center Relevant','Quadrant Relevant'};
        rel_save_name = 'cent_quad';

    end

%Irrelevant stimuli
elseif isequal(relevant_condition, 'irrel')
    
    %Condition name
    condition_name = 'irrelevant';
    
    %Center Irrelevant
    if isequal(location_set, 'c')

        location_cell = {'Center Irrelevant'};
        rel_save_name = 'cent';

    %Quadrant Irrelevant
    elseif isequal(location_set, 'q')

        location_cell = {'Quadrant Irrelevant'};
        rel_save_name = 'quad';

    %Combine Center and Quadrant Irelevant   
    elseif isequal(location_set, 'combine')

        location_cell = {'Center Irrelevant','Quadrant Irrelevant'};
        rel_save_name = 'cent_quad';

    end
    
% Combine Relevant and Irrelevant stimuli
elseif isequal(relevant_condition, 'combine')
    
    %Condition name
    condition_name = 'rel_irrel';
    
    %Center Irrelevant
    if isequal(location_set, 'c')

        location_cell = {'Center Relevant', 'Center Irrelevant'};
        rel_save_name = 'cent';

    %Quadrant Irrelevant
    elseif isequal(location_set, 'q')

        location_cell = {'Quadrant Relevant', 'Quadrant Irrelevant'};
        rel_save_name = 'quad';

    %Both Center and Quadrant Irelevant   
    elseif isequal(location_set, 'combine')

        location_cell = {'Center Relevant','Quadrant Relevant','Center Irrelevant','Quadrant Irrelevant'};
        rel_save_name = 'cent_quad';

    end    

end

%Imaging modality
prompt_4 = 'Modality type [EEG, MRI, both]: ';
modality = input(prompt_4,'s');

%EEG
if isequal(modality,'EEG')
   
    modality_name = 'EEG';
    modality_cell = {'EEG'};
    
%MRI
elseif isequal(modality,'MRI')
    
    modality_name = 'MRI';
    modality_cell = {'MRI'};
    
%EEG and MRI
elseif isequal(modality,'both')
    
    modality_name = 'EEG_MRI';
    modality_cell = {'EEG','MRI'};
    
end

%Remove trials with blink at stimulus onset 
prompt_5 = 'Remove trials with blink at stim [y,n]: ';
remove_blink_trials = input(prompt_5,'s');

%% Local directories
if isequal(run_location, 's')
    
    %Add group directory to path
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Data/EEG')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/Behavioral Analysis')
    addpath(genpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/No Report Paradigm Machine Learning/Training Testing and Plotting Functions'))
    
    %Subject EEG and MRI dir
    EEG_sub_dir = '/mnt/Data28/HNCT No Report Paradigm/Subject Analysis EEG';
    MRI_sub_dir = '/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI';
    
    %Save directory
    save_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Data/EyeLink';
        
elseif isequal(run_location, 'l')
       
end

%% Aggregate Eyelink Data Over Subjects

%Initialized group variables
group_pupil_data = [];
group_blink_data = [];
group_microsac_data = [];
group_subject_trial_index = [];
group_stacked_SVM_avg_scores = [];

%Loop over stimulus locations
for loc = 1:length(location_cell)

    %Select location
    current_location = location_cell{loc};

    %Loop over modality type
    for mod = 1:length(modality_cell)

        %Select current modality
        current_modality = modality_cell{mod};

        %Subject eyelink data path 
        if isequal(run_location, 'l')

        elseif isequal(run_location, 's')

            %EEG subject directories
            if isequal(current_modality, 'EEG')
                
                %Threshold
                if isequal(opacity,'threshold')
                    
                    %Model path
                    model_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/Machine Learning Analysis/EEG Session'...
                        ,current_location,'Base_boosted_tree_lin_gaus_SVM_stacked_lin_SVM');
                
                %Blank
                elseif isequal(opacity,'blank')
                    
                    %Model path
                    model_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/Machine Learning Analysis/EEG Session'...
                        ,current_location,'Base_boosted_tree_lin_gaus_SVM_stacked_lin_SVM_blank_stim');
                                        
                end
                
            %MRI subject directories
            elseif isequal(current_modality, 'MRI')
                
                %Threshold
                if isequal(opacity,'threshold')
                    
                    %Model path
                    model_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/Machine Learning Analysis/MRI Session'...
                        ,current_location,'Base_boosted_tree_lin_gaus_SVM_stacked_lin_SVM');
                
                %Blank
                elseif isequal(opacity,'blank')
                    
                    %Model path
                    model_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/Machine Learning Analysis/MRI Session'...
                        ,current_location,'Base_boosted_tree_lin_gaus_SVM_stacked_lin_SVM_blank_stim');

                end
                
            end

        end

        %% Store Model Confidence Scores

        %Load model direcotry 
        cd(model_dir)
        
        %Relevant data model with threshold opacity
        %if isequal(relevant_condition, 'rel')
        if isequal(opacity,'threshold') && isequal(relevant_condition,'rel')
            
            load('Classifier_models_scores_labels.mat')
            
            %Store the avg confidence scores into group variable
            group_stacked_SVM_avg_scores = cat(1,group_stacked_SVM_avg_scores,stacked_scores);
        
        %Relevant and blank; Irrelevant threshold and blank data model
        else %isequal(relevant_condition, 'irrel')
            
            %Unique file name for blank stimuli
            if isequal(opacity,'blank')
                
                load('Classifier_scores_labels_blank_stim.mat')
                
            else
                
                load('Classifier_noreport_scores_labels.mat')
                
            end
            
            %Store the avg confidence scores into group variable
            group_stacked_SVM_avg_scores = cat(1,group_stacked_SVM_avg_scores,stacked_SVM_avg_scores);
        
        end

        %% Store Eyelink Data

        %Load subject EEG and EEG-eyelink data and add to a group variable. Output
        %the group extracted eyelink, EEG, and true CP/CnP labels. 
        if isequal(current_modality, 'EEG') && ismember(current_location,{'Quadrant Irrelevant','Center Irrelevant'}) %&& isequal(relevant_condition,'irrel')

            disp('Extract Irrelevant EEG Eyelink Data')

            %Extract the irrelevant eyelink trials - EEG
            [noreport_data_pupil, noreport_data_tsfreshpupil, noreport_data_Ygaze, noreport_data_Xgaze, noreport_data_blink, noreport_data_microsac,...
                group_subject, stimulus_opacity, subject_trial_index] = combine_eyelink_EEG_irrelevant_epochs_function(EEG_sub_dir,current_location, run_location);

        elseif isequal(current_modality, 'MRI') && ismember(current_location,{'Quadrant Irrelevant','Center Irrelevant'}) %&& isequal(relevant_condition,'irrel')

            disp('Extract Irrelevant MRI Eyelink Data')

            %Extract the relevant eyelink trials - MRI
            [noreport_data_pupil, noreport_data_tsfreshpupil, noreport_data_Ygaze, noreport_data_Xgaze, noreport_data_blink, noreport_data_microsac, group_subject,...
                group_run_number, group_trial_number, group_condition, group_relevant_location, group_BOLD_epoch_count, number_of_MRI_subjects, subject_trial_index] = ...
                combine_eyelink_MRI_irrelevant_epochs_function(MRI_sub_dir, current_location, run_location, opacity);

        elseif isequal(current_modality, 'EEG') && ismember(current_location,{'Quadrant Relevant','Center Relevant'}) %&& isequal(relevant_condition,'rel')

            disp('Extract Relevant EEG Eyelink Data')

            %Extract the relevant eyelink trials - EEG
            [report_data_pupil, report_data_tsfreshpupil, report_data_Ygaze, report_data_Xgaze, report_data_blink, report_data_microsac, ...
                true_labels, stimulus_opacity, group_subject, group_EEG_epoch_count, number_of_EEG_subjects, subject_trial_index] = ...
                combine_eyelink_EEG_relevant_epochs_function(EEG_sub_dir, current_location, run_location);

        elseif isequal(current_modality, 'MRI') && ismember(current_location,{'Quadrant Relevant','Center Relevant'}) %&& isequal(relevant_condition,'rel')

            disp('Extract Relevant MRI Eyelink Data')

            %Extract the relevant eyelink trials  - MRI
            [report_data_pupil, report_data_tsfreshpupil, report_data_Ygaze, report_data_Xgaze, report_data_blink, report_data_microsac, true_labels, ...
                group_subject, group_run_number, group_trial_number, group_condition, group_relevant_location, group_BOLD_epoch_count, number_of_MRI_subjects, subject_trial_index] = ...
                combine_eyelink_MRI_relevant_epochs_function(MRI_sub_dir, current_location, run_location, opacity);

        end

        %Store group eyelink data
        if ismember(current_location,{'Quadrant Relevant','Center Relevant'})

            group_pupil_data = cat(1,group_pupil_data,report_data_pupil);
            group_blink_data = cat(1,group_blink_data,report_data_blink);
            group_microsac_data = cat(1,group_microsac_data,report_data_microsac);

        elseif ismember(current_location,{'Quadrant Irrelevant','Center Irrelevant'})

            group_pupil_data = cat(1,group_pupil_data,noreport_data_pupil);
            group_blink_data = cat(1,group_blink_data,noreport_data_blink);
            group_microsac_data = cat(1,group_microsac_data,noreport_data_microsac);

        end

        %Store subject trial index
        group_subject_trial_index = cat(1,group_subject_trial_index,subject_trial_index);
        
        %Clear variables
        clearvars report* noreport* subject_trial_index

    end

end

%% Confirm correct number of trials

if not(isequal(length(group_subject_trial_index),length(group_stacked_SVM_avg_scores)))
    
   error('Mismatched group trial number!')
    
end
    
%% Find PP and PnP Trials and Create Subject Matrix

%Confidence Score Thresholds
confidence_score_thresholds = [0:0.25:1];   

%Loop over confidence scores
for thres = 1:length(confidence_score_thresholds)

    %Reset group variables
    group_PP_pupil_data = [];
    group_PnP_pupil_data = [];
    group_PP_blink_data = [];
    group_PnP_blink_data = [];
    group_PP_microsac_data = [];
    group_PnP_microsac_data = [];

    %Define threshold 
    current_threshold = confidence_score_thresholds(thres);

    disp(['Finding PP and PnP trials for score threshold ',num2str(current_threshold)])

    %Extract the predicted class trials
    PP_pupil = group_pupil_data(group_stacked_SVM_avg_scores(:,1) > current_threshold,:);
    PnP_pupil = group_pupil_data(group_stacked_SVM_avg_scores(:,1) < -current_threshold,:);

    PP_blink = group_blink_data(group_stacked_SVM_avg_scores(:,1) > current_threshold,:);
    PnP_blink = group_blink_data(group_stacked_SVM_avg_scores(:,1) < -current_threshold,:);

    PP_microsac = group_microsac_data(group_stacked_SVM_avg_scores(:,1) > current_threshold,:);
    PnP_microsac = group_microsac_data(group_stacked_SVM_avg_scores(:,1) < -current_threshold,:);

    PP_subject_index = group_subject_trial_index(group_stacked_SVM_avg_scores(:,1) > current_threshold,:);
    PnP_subject_index = group_subject_trial_index(group_stacked_SVM_avg_scores(:,1) < -current_threshold,:);

    PP_SVM_scores = group_stacked_SVM_avg_scores(group_stacked_SVM_avg_scores(:,1) > current_threshold,:);
    PnP_SVM_scores = group_stacked_SVM_avg_scores(group_stacked_SVM_avg_scores(:,1) < -current_threshold,:);

    %Remove trials with blinks at stim onset
    if isequal(remove_blink_trials,'y')

        %Find trial with blink at stimulus presenation
        PP_blink_at_stim_trials = find(sum(PP_blink(:,6001:6050),2) > 0);
        PnP_blink_at_stim_trials = find(sum(PnP_blink(:,6001:6050),2) > 0);
        
        %Remove trials from Eyelink, subject_index, and score matrices
        PP_pupil(PP_blink_at_stim_trials,:) = [];
        PP_blink(PP_blink_at_stim_trials,:) = [];
        PP_microsac(PP_blink_at_stim_trials,:) = [];
        PP_subject_index(PP_blink_at_stim_trials,:) = [];
        PP_SVM_scores(PP_blink_at_stim_trials,:) = [];

        PnP_pupil(PnP_blink_at_stim_trials,:) = [];
        PnP_blink(PnP_blink_at_stim_trials,:) = [];
        PnP_microsac(PnP_blink_at_stim_trials,:) = [];
        PnP_subject_index(PnP_blink_at_stim_trials,:) = [];
        PnP_SVM_scores(PnP_blink_at_stim_trials,:) = [];

        %Count the number of blink at stim trials
        PP_blink_rejected_trial_num = length(PP_blink_at_stim_trials);
        PnP_blink_rejected_trial_num = length(PnP_blink_at_stim_trials);

    end
    
    %Count the number of PP and PnP trials
    PP_total_trial_num = sum(group_stacked_SVM_avg_scores(:,1) > current_threshold);
    PnP_total_trial_num = sum(group_stacked_SVM_avg_scores(:,1) < -current_threshold);

    %Average trials within subjects

    %Unique subject IDs
    PP_epochs_subjects_list = unique(PP_subject_index);
    PnP_epochs_subjects_list = unique(PnP_subject_index);

    %Loop over PP subjects 
    for sub = 1:length(PP_epochs_subjects_list)

        %Initialize variable
        subject_PP_pupil_data = [];
        subject_PnP_pupil_data = [];
        subject_PP_blink_data = [];
        subject_PnP_blink_data = [];
        subject_PP_microsac_data = [];
        subject_PnP_microsac_data = [];

        %Define current ID in loop
        current_ID = PP_epochs_subjects_list{sub};

        %Find all trails belonging to current ID
        current_ID_PP_trial_index = find(strcmp(PP_subject_index,current_ID));

        %Average PP for current ID
        subject_PP_pupil_data = nanmean(PP_pupil(current_ID_PP_trial_index,:),1);
        subject_PP_blink_data = nanmean(PP_blink(current_ID_PP_trial_index,:),1);
        subject_PP_microsac_data = nanmean(PP_microsac(current_ID_PP_trial_index,:),1);

        % Create group matrix of [time x subjects] - PP
        if isempty(group_PP_pupil_data) 

            group_PP_pupil_data(:,1) = subject_PP_pupil_data; 
            group_PP_blink_data(:,1) = subject_PP_blink_data;
            group_PP_microsac_data(:,1) = subject_PP_microsac_data;

        else

            group_PP_pupil_data(:,(end+1)) = subject_PP_pupil_data; 
            group_PP_blink_data(:,(end+1)) = subject_PP_blink_data;
            group_PP_microsac_data(:,(end+1)) = subject_PP_microsac_data;

        end 

    end

    %Loop over PnP subjects 
    for sub = 1:length(PnP_epochs_subjects_list)

        %Define current ID in loop
        current_ID = PnP_epochs_subjects_list{sub};

        %Find all trails belonging to current ID
        current_ID_PnP_trial_index = find(strcmp(PnP_subject_index,current_ID));

        %Average PnP for current ID
        subject_PnP_pupil_data = nanmean(PnP_pupil(current_ID_PnP_trial_index,:),1);
        subject_PnP_blink_data = nanmean(PnP_blink(current_ID_PnP_trial_index,:),1);
        subject_PnP_microsac_data = nanmean(PnP_microsac(current_ID_PnP_trial_index,:),1);

        % Create group matrix of [time x subjects] - PnP
        if isempty(group_PnP_pupil_data) 

            group_PnP_pupil_data(:,1) = subject_PnP_pupil_data; 
            group_PnP_blink_data(:,1) = subject_PnP_blink_data;
            group_PnP_microsac_data(:,1) = subject_PnP_microsac_data;

        else

            group_PnP_pupil_data(:,(end+1)) = subject_PnP_pupil_data; 
            group_PnP_blink_data(:,(end+1)) = subject_PnP_blink_data;
            group_PnP_microsac_data(:,(end+1)) = subject_PnP_microsac_data;

        end 

    end  

    %Save group EyeLink data
    cd(save_dir)
    
    if isequal(remove_blink_trials,'y')
        
        save(['Group_eyelink_',condition_name,'_PP_PnP_',rel_save_name,'_',opacity,'_',modality_name,'_score_thres_',num2str(current_threshold),'_data.mat'],'PP_epochs_subjects_list','PnP_epochs_subjects_list', ...
            'group_PP_pupil_data','group_PnP_pupil_data','group_PP_blink_data','group_PnP_blink_data','group_PP_microsac_data','group_PnP_microsac_data',...
            'PP_total_trial_num','PnP_total_trial_num','PP_blink_rejected_trial_num','PnP_blink_rejected_trial_num','-v7.3')

    else
        
        save(['Group_eyelink_',condition_name,'_PP_PnP_',rel_save_name,'_',opacity,'_',modality_name,'_score_thres_',num2str(current_threshold),'_no_blink_rejections_data.mat'],'PP_epochs_subjects_list','PnP_epochs_subjects_list', ...
            'group_PP_pupil_data','group_PnP_pupil_data','group_PP_blink_data','group_PnP_blink_data','group_PP_microsac_data','group_PnP_microsac_data',...
            'PP_total_trial_num','PnP_total_trial_num','-v7.3')
        
    end
    
end
