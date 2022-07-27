%% Combine PC BOLD Data Across Paradigm Dataset

%This code will load the report and noreport BOLD PC data for CP and CnP
%trials and combine over the subject dimension and save in the combined
%data directory.

%Written: Sharif Kronemer
%Date: 5/18/2021

clear

%% Directories 

%Report dataset group data directory
noreport_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Data/MRI/BOLD PC';

%No report dataset group data directory
report_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Group Data/BOLD PC';

%Save dir
save_dir = '/mnt/Data8/HNCT NRP and Report Only Paradigm Analyses/Group Data/MRI/BOLD PC';

%% Combine report and no report datasets

%Report data variable
noreport_file = 'Group_PC_BOLD_CP_CnP_cent_quad_thres_data.mat';

%No report data variable
report_file = 'Group_PC_BOLD_CP_CnP_poststim_15s_data.mat';

%Load report data
load(fullfile(report_dir,report_file))

report_group_CP_BOLD_PC_data = group_CP_BOLD_PC_data;
report_group_CnP_BOLD_PC_data = group_CnP_BOLD_PC_data;
report_all_included_subjects_list = all_included_subjects_list;

clearvars group* 

%Load noreport data
load(fullfile(noreport_dir,noreport_file))

noreport_group_CP_BOLD_PC_data = group_CP_BOLD_PC_data;
noreport_group_CnP_BOLD_PC_data = group_CnP_BOLD_PC_data;
noreport_all_included_subjects_list = all_included_subjects_list;

clearvars group* 

%Combine datasets over the subject dimension
group_CP_BOLD_PC_data = cat(3, report_group_CP_BOLD_PC_data, noreport_group_CP_BOLD_PC_data);

clearvars report_group_CP* noreport_group_CP*

group_CnP_BOLD_PC_data = cat(3, report_group_CnP_BOLD_PC_data, noreport_group_CnP_BOLD_PC_data);

clearvars report_group_CnP* noreport_group_CnP*

all_included_subjects_list = [report_all_included_subjects_list;noreport_all_included_subjects_list];

%% Save combined data

%Save combined data
cd(save_dir)
save('Group_PC_BOLD_CP_CnP_threshold_data.mat','group_CP_BOLD_PC_data','group_CnP_BOLD_PC_data',...
    'all_included_subjects_list','-v7.3')
