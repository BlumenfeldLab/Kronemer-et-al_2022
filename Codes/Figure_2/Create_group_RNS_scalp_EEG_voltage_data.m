%% Create Group 10-20 Scalp EEG RNS 

%This code aggregates scalp EEG epochs across RNS subjects into a group
%matrix of channels x time x subjects

%Written by: Sharif Kronemer
%Date: 5/26/2021

clear

%% Directories 

%Reference type
reference_type = 'Average Reference';

%Group data save dir
group_dir = fullfile('R:\RNS Study\Analysis\Group Data',reference_type);
mkdir(group_dir)

% Subject List
subject_list = {'471MH','489RD','490KB','491GS','535BP','536JD'};

%% Aggregate across subjects

%Group variables - good trials (not rejected in preprocessing)
group_CP_good_EEG_voltage_data = [];
group_CnP_good_EEG_voltage_data = [];

%Group variables - all trials 
group_CP_all_EEG_voltage_data = [];
group_CnP_all_EEG_voltage_data = [];

%Total trial count
CP_good_trial_num = 0;
CnP_good_trial_num = 0;
CP_all_trial_num = 0;
CnP_all_trial_num = 0;

%Loop over subjects
for sub = 1:length(subject_list)

    %Select subject ID
    ID = subject_list{sub};
    
    disp(['Running subject ', ID])
    
    %Data matrix
    data_dir = fullfile('R:\RNS Study\Analysis\Subject Analysis\',ID,'\Scalp EEG Analysis\Preprocessed Data',reference_type,'Voltage Epochs');
    
    %Load CP and CnP
    load(fullfile(data_dir,'CP_4s_epochs.mat'))
    load(fullfile(data_dir,'CnP_4s_epochs.mat'))
    
    %Add subject data to group matrix [channel x time x subjects]
    group_CP_good_EEG_voltage_data = cat(3,group_CP_good_EEG_voltage_data,nanmean(CP_good_epochs,3));
    group_CnP_good_EEG_voltage_data = cat(3,group_CnP_good_EEG_voltage_data,nanmean(CnP_good_epochs,3));
    group_CP_all_EEG_voltage_data = cat(3,group_CP_all_EEG_voltage_data,nanmean(CP_all_epochs,3));
    group_CnP_all_EEG_voltage_data = cat(3,group_CnP_all_EEG_voltage_data,nanmean(CnP_all_epochs,3));
    
    %Add to trial count
    CP_good_trial_num = CP_good_trial_num + size(CP_good_epochs,3);
    CnP_good_trial_num = CnP_good_trial_num + size(CnP_good_epochs,3);
    CP_all_trial_num = CP_all_trial_num + size(CP_all_epochs,3);
    CnP_all_trial_num = CnP_all_trial_num + size(CnP_all_epochs,3);
    
end

%Save group data
cd(group_dir)
save('Group_voltage_scalp_EEG.mat', 'group_CP*','group_CnP*','CP_all_trial_num','CnP_all_trial_num','CP_good_trial_num','CnP_good_trial_num')
