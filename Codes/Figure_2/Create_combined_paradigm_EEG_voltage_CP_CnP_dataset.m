%% Combine EEG data across report and no-report datasets

%Note: You must manually select the filenames to combine.

%Written by: Sharif Kronemer
%Date: 5/19/2021

clear

%% Directories

%Report directory
report_dir = '/mnt/Data6/HNCT sEEG Study/sEEG Analysis/Group Data/Voltage';

%No report directory
noreport_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Data/EEG/Voltage';

%Save directory
save_dir = '/mnt/Data8/HNCT NRP and Report Only Paradigm Analyses/Group Data/EEG/Voltage';

%Define the files to load 
report_filename = 'Group_EEG_voltage_CP_CnP_1_15s_data.mat'; 
noreport_filename = 'Group_EEG_voltage_CP_CnP_quad_threshold_data.mat';
save_filename = 'Group_EEG_voltage_CP_CnP_quad_threshold_1_15s_data.mat';

%% Load EEG data

%Load report data
load(fullfile(report_dir,report_filename))

%Rename variables
report_CP = group_CP_EEG_voltage_data;
report_CnP = group_CnP_EEG_voltage_data;
report_sub_list_CP = CP_epochs_subjects_list;
report_sub_list_CnP = CnP_epochs_subjects_list;

clearvars group* CP_epochs_subjects_list CnP_epochs_subjects_list

%Load noreport data
load(fullfile(noreport_dir,noreport_filename))

%Rename variables
noreport_CP = group_CP_EEG_voltage_data;
noreport_CnP = group_CnP_EEG_voltage_data;
noreport_sub_list_CP = CP_epochs_subjects_list;
noreport_sub_list_CnP = CnP_epochs_subjects_list;

clearvars group* CP_epochs_subjects_list CnP_epochs_subjects_list

%% Combine CP and CnP data

disp('***Combining data***')

%CP and CnP EEG voltage 
group_CP_EEG_voltage_data = cat(3,report_CP,noreport_CP);
group_CnP_EEG_voltage_data = cat(3,report_CnP,noreport_CnP);

%Subject list
CP_epochs_subjects_list = [report_sub_list_CP; noreport_sub_list_CP];
CnP_epochs_subjects_list = [report_sub_list_CnP; noreport_sub_list_CnP];

%% Save Data

cd(save_dir)
save(save_filename,'group_CP_EEG_voltage_data','group_CnP_EEG_voltage_data',...
    'CP_epochs_subjects_list','CnP_epochs_subjects_list','-v7.3')
