%% All-in-one Percent Change Analysis 

%This script is design to take preprocessed nifiti files and extract BOLD
%signal within designated events. 

%(1) Finds the CP and CnP and irrelevant events from log file
%(2) Extracts BOLD signal 
%(3) Artifact rejection and percent change BOLD signal calculation
%(4) Percent Change Epoch Cutting, Binning, and Normalization

%Written by: Sharif I. Kronemer
%Date: 10/14/2020
%Modified: 5/8/2021

clear

%Note: The script could be updated to run a single specified subject ID.
%The current version of the script finds all subject IDs and runs all
%subjects. 

% %% Select Subject Data to Analyze
% prompt1 = 'Subject ID: '; %prompt ID
% ID = input(prompt1, 's');

%% Location Relevant

prompt2 = 'Relevant location (q/c/both): '; %prompt ID
relevant = input(prompt2, 's');

%% Perception Response

prompt3 = 'Perception response yes (1/2): '; 
Yes = input(prompt3);

%Define relevant directory
if isequal(relevant, 'q')
    
   relevant_cell = {'Quadrant Relevant'};
   irrelevant_cell = {'Center Irrelevant'};
       
elseif isequal(relevant, 'c')
    
   relevant_cell = {'Center Relevant'};
   irrelevant_cell = {'Quadrant Irrelevant'};
   
elseif isequal(relevant, 'both')
    
   relevant_cell = {'Quadrant Relevant','Center Relevant'};
   irrelevant_cell = {'Center Irrelevant', 'Quadrant Irrelevant'};
       
end

%% Run location

prompt4 = 'Run location (l/s): '; 
run_location = input(prompt4, 's');

%Subject list
if isequal(run_location, 's')
    
    subject_list = dir('/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI');

elseif isequal(run_location, 'l')
    
    subject_list = dir('Y:\HNCT No Report Paradigm\Subject Analysis MRI');

end

%Loop over subjects
for sub = 1:length(subject_list) 

    %Define ID
    ID = subject_list(sub).name;

    %% Cycle Over Relevant Locations

    tic

    %Loop over relevant locations
    for rel = 1:size(relevant_cell,2)

        %Select the irrelevant_dir and relevant_dir
        relevant_dir = relevant_cell{rel};
        irrelevant_dir = irrelevant_cell{rel};

        %Specific subject rejection - all runs rejected    
        if ismember(ID,{'567','568','648'}) && isequal(relevant_dir,'Center Relevant')
                       
            disp(['All runs rejected by motion - Skipping ',ID])
            continue
            
        %Specific subject rejection - all runs rejected    
        elseif ismember(ID,{'632'}) && isequal(relevant_dir,'Quadrant Relevant')
                       
            disp(['All runs rejected by motion - Skipping ',ID])
            continue  
            
        else
            
            disp(['Running - ',ID, ' ',relevant_dir, ' and ',irrelevant_dir])

        end
        
        %% Global Directories and Paths

        %Local
        if isequal(run_location, 'l')

            %rootpath = fullfile('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Subject Analysis\',ID,'Perception Task',relevant_dir,'\MRI Session\MRI Analysis\Extracted voxel data');
            rootpath = fullfile('Y:\HNCT No Report Paradigm\Subject Analysis MRI\',ID,'Perception Task',relevant_dir,'\MRI Session\MRI Analysis\Extracted voxel data');
            
            irrel_rootpath = fullfile('Y:\HNCT No Report Paradigm\Subject Analysis MRI',ID,'Perception Task',irrelevant_dir,'MRI Analysis\Extracted voxel data');
            
            %beh_data_path = fullfile('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Subject Analysis\',ID,'Perception Task',relevant_dir,'\MRI Session\MRI Analysis\Event Times');
            beh_data_path = fullfile('Y:\HNCT No Report Paradigm\Subject Analysis MRI\',ID,'Perception Task',relevant_dir,'\MRI Session\MRI Analysis\Event Times');

            %subject_dir = 'S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Subject Analysis';
            subject_dir = 'Y:\HNCT No Report Paradigm\Subject Analysis MRI';

            %Add analysis folders to path
            addpath(genpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\Raw Timecourses\Supplementary functions'))
            addpath(genpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\Extract voxel data\Supplementary functions'))
            addpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\spm12_radiological')
            addpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\MR Behavioral Analysis')
            addpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\Extract voxel data')
            addpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\marsbar-0.44')
            addpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\Percent Change Analysis\Subject Level')

            %Behavior data and save folders
            beh_folder = fullfile('M:\Subject Raw Data\', ID, 'Perception Task', relevant_dir,'MRI Session\Behavioral Data');

            %savefolder = fullfile('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Subject Analysis\',ID, 'Perception Task', relevant_dir,'\MRI Session\MRI Analysis\Event Times');
            savefolder = fullfile('Y:\HNCT No Report Paradigm\Subject Analysis MRI\',ID, 'Perception Task', relevant_dir,'\MRI Session\MRI Analysis\Event Times');

            %Motion artifact directory
            motion_dir = 'S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Group Analysis\MRI Analysis\Movement Artifact';
            
            irrel_savefolder = fullfile('Y:\HNCT No Report Paradigm\Subject Analysis MRI',ID, 'Perception Task', irrelevant_dir,'MRI Analysis\Event Times');

        %Server
        elseif isequal(run_location, 's')

            %rootpath = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Subject Analysis',ID,'Perception Task',relevant_dir,'MRI Session/MRI Analysis/Extracted voxel data');
            rootpath = fullfile('/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI',ID,'Perception Task',relevant_dir,'MRI Session/MRI Analysis/Extracted voxel data');

            irrel_rootpath = fullfile('/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI',ID,'Perception Task',irrelevant_dir,'MRI Analysis/Extracted voxel data');

            %beh_data_path = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Subject Analysis',ID,'Perception Task',relevant_dir,'MRI Session/MRI Analysis/Event Times');
            beh_data_path = fullfile('/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI',ID,'Perception Task',relevant_dir,'MRI Session/MRI Analysis/Event Times');
            irrel_beh_data_path = fullfile('/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI',ID,'Perception Task',irrelevant_dir,'MRI Analysis/Event Times');

            %subject_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Subject Analysis';
            subject_dir = '/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI';

            %Add analysis folders to path
            addpath(genpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Raw Timecourses/Supplementary functions'))
            addpath(genpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Extract voxel data/Supplementary functions'))
            addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/spm12_radiological')
            addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/MR Behavioral Analysis')
            addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Extract voxel data')
            addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/marsbar-0.44')
            addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Percent Change Analysis/Subject Level')

            %Behavior data and save folders
            beh_folder = fullfile('/mnt/Data15/Subject Raw Data/', ID, 'Perception Task', relevant_dir,'MRI Session/Behavioral Data');

            %savefolder = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Subject Analysis/',ID, 'Perception Task', relevant_dir,'/MRI Session/MRI Analysis/Event Times');
            savefolder = fullfile('/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI/',ID, 'Perception Task', relevant_dir,'/MRI Session/MRI Analysis/Event Times');

            irrel_savefolder = fullfile('/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI/',ID, 'Perception Task', irrelevant_dir,'MRI Analysis/Event Times');

            %Motion artifact directory
            motion_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/MRI Analysis/Movement Artifact';

        end

        %Preprocessing nii files path
        preprocessing_dir = fullfile(subject_dir, ID, 'Perception Task', relevant_dir, 'MRI Session', 'MRI Analysis', 'Preprocessed Images');

        %If MR data does not exist - Skip
        if not(exist(preprocessing_dir))

            disp(['Data does not exist - Skipping ',ID, ' ', relevant_dir ' and ', irrelevant_dir])
            continue

        end

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
        %Find_CP_CnP_Events_From_Log_File_Function(Yes, beh_folder, savefolder);
        
        %Blank Stim - Relevant Condition
        %Find_Blank_Events_From_Log_File_Function(beh_folder, savefolder);

        %Irrelevant Stimuli Events
        %Inputs: beh_folder and savefolder
        %Find_Task_Irrelevant_Events_From_Log_File_Function(beh_folder, savefolder)
        
        %Jitter Events
        Find_Jitter_Events_From_Log_File_Function(beh_folder, savefolder, ID, relevant_dir)

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

        % Stimulus centered epochs
        %Percent_change_BOLD_epoch_cutting_function(ID, rootpath, irrel_rootpath, run_num, run_list, relevant_dir, motion_dir, savefolder, irrel_savefolder, beh_data_path, run_location);
              
        % Jitter centered epochs
        %OLD METHOD: Percent_change_BOLD_epoch_cutting_shifted_function(ID, rootpath, irrel_rootpath, run_num, run_list, relevant_dir, irrelevant_dir, motion_dir, savefolder, irrel_savefolder, beh_data_path, irrel_beh_data_path, run_location);
        Percent_change_BOLD_jitter_epoch_cutting_function(ID, rootpath, run_num, run_list, relevant_dir, motion_dir, savefolder, beh_data_path, run_location);

    end

    toc

end