%% All-in-one Percent Change Analysis - Report only paradigm

%***NOTE***: Not all elements of this code has been updated for the report
%only paradigm, which is adapted from the no report paradigm version of
%this code.

%This script is design to take preprocessed nifiti files and extract BOLD
%signal within designated events. 

%(1) Finds the CP and CnP and irrelevant events from log file
%(2) Extracts BOLD signal 
%(3) Artifact rejection and percent change BOLD signal calculation
%(4) Percent Change Epoch Cutting, Binning, and Normalization

%Written by: Sharif I. Kronemer
%Date: 10/14/2020
%Modified: 4/26/2021

clear

%% Select Subject Data to Analyze
prompt1 = 'Subject ID: '; %prompt ID
ID = input(prompt1, 's');

%% Run location

prompt4 = 'Run location (l/s): '; 
run_location = input(prompt4, 's');

%% Global Directories and Paths

%Local
if isequal(run_location, 'l')

    rootpath = fullfile('Z:\Sharif Data F\HNCT fMRI Study\fMRI Analysis\Subject Analysis\',ID,'\MRI Analysis\Extracted voxel data');

    beh_data_path = fullfile('Z:\Sharif Data F\HNCT fMRI Study\fMRI Analysis\Subject Analysis\',ID,'\Behavioral Analysis');

    subject_dir = 'Z:\Sharif Data F\HNCT fMRI Study\fMRI Analysis\Subject Analysis';

    %Add analysis folders to path
    addpath(genpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\Raw Timecourses\Supplementary functions'))
    addpath(genpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\Extract voxel data\Supplementary functions'))
    addpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\spm12_radiological')
    addpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\MR Behavioral Analysis')
    addpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\Extract voxel data')
    addpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\marsbar-0.44')
    addpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\Percent Change Analysis\Subject Level')

    %Behavior data and save folders
    %beh_folder = fullfile('M:\Subject Raw Data\', ID,'Behavioral Data');

    %savefolder = fullfile('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Subject Analysis\',ID, 'Perception Task', relevant_dir,'\MRI Session\MRI Analysis\Event Times');
    %savefolder = fullfile('Y:\HNCT No Report Paradigm\Subject Analysis MRI\',ID, 'Perception Task', relevant_dir,'\MRI Session\MRI Analysis\Event Times');

    %Motion artifact directory
    motion_dir = 'Z:\Sharif Data F\HNCT fMRI Study\fMRI Analysis\Group Analysis\Movement Artifact';

%Server
elseif isequal(run_location, 's')

    rootpath = fullfile('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Subject Analysis/',ID,'/MRI Analysis/Extracted voxel data');

    beh_data_path = fullfile('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Subject Analysis/',ID,'/Behavioral Analysis');

    subject_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Subject Analysis';

    %Add analysis folders to path
    addpath(genpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Raw Timecourses/Supplementary functions'))
    addpath(genpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Extract voxel data/Supplementary functions'))
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/spm12_radiological')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/MR Behavioral Analysis')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Extract voxel data')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/marsbar-0.44')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Percent Change Analysis/Subject Level')

    addpath('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Percent Change Maps/Percent Change Maps')
    
    %Behavior data and save folders
    %beh_folder = fullfile('/mnt/Data15/Subject Raw Data/', ID, 'Perception Task', relevant_dir,'MRI Session/Behavioral Data');

    %savefolder = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Subject Analysis/',ID, 'Perception Task', relevant_dir,'/MRI Session/MRI Analysis/Event Times');
    %savefolder = fullfile('/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI/',ID, 'Perception Task', relevant_dir,'/MRI Session/MRI Analysis/Event Times');

    %Motion artifact directory
    %motion_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/MRI Analysis/Movement Artifact';
    motion_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Group Analysis/Movement Artifact';

end

%Preprocessing nii files path
preprocessing_dir = fullfile(subject_dir,ID,'MRI Analysis','Preprocessing');

%Create directory of runs
cd(preprocessing_dir)
run_list = dir('Run*');
run_num = length(run_list); %Find the number of runs completed

%% STEP 1 - Mine Log File for Relevant and Irrelevant Events

disp('**STEP 1 - Extract Events')

%Make the savefolder if it does not exist
if ~exist(savefolder)

    disp('Making save directory for behavioral event times')
    mkdir(savefolder)

end

%CP and CnP Events
%Inputs: Yes (perception button press), beh_folder, and savefolder
Find_CP_CnP_Events_From_Log_File_Function(Yes, beh_folder, savefolder);

%Irrelevant Stimuli Events
%Inputs: beh_folder and savefolder
Find_Task_Irrelevant_Events_From_Log_File_Function(beh_folder, savefolder)

%% STEP 2 - Extract BOLD Signal

disp('**STEP 2 - Extract BOLD Signal')

%Loop over runs
for run = 1:run_num

    %disp(['Running - Run ',num2str(run)])
    disp(['Running - ',run_list(run).name])


    % Define the data directory
    dataDir = fullfile(preprocessing_dir,run_list(run).name);

    % Mask out gray-matter
    if isequal(run_location, 's')

        roiPath = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/marsbar-aal-0.2/Gray Matter with hypo NB SN subthal';

    elseif isequal(run_location, 'l')

        roiPath = 'Z:\Sharif Data F\HNCT fMRI Study\fMRI Analysis\Analysis Codes\marsbar-aal-0.2\Gray Matter with hypo NB SN subthal';

    end

    % Load gray matter mask
    roiFile = 'Gray_Matter_with_hypo_NB_SN_subthal_roi.mat';
    roiPath = fullfile(roiPath, roiFile);

    % Define results directory
    resultsDir = fullfile(subject_dir,ID,'Perception Task', relevant_dir,'MRI Session','MRI Analysis','Extracted voxel data',run_list(run).name); 

    % Make results directory
    if ~exist(resultsDir,'dir')

        mkdir(resultsDir)

    end

    % Run BOLD extraction
    extractfMRIdata_fun(resultsDir, dataDir, roiPath, run_location);

end

%% STEP 3 - Artifact Rejection and Calculate Percent Change

disp('**STEP 3 - Artifact Rejection')

Percent_Change_Map_Code_jitter_norm_new_bin_filter_Function(run_location, relevant_dir, ID, rootpath, run_num, run_list)

%% STEP 4 - Percent Change Epoch Cutting, Binning, and Normalization

disp('**STEP 4 - Cutting, Binning, and Normalizing Percent Change Epochs')

Percent_change_BOLD_epoch_cutting_report_only_paradigm_function(ID, rootpath, beh_data_path, run_num, run_list, motion_dir);
