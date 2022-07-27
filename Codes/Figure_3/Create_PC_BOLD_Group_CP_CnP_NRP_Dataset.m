%% Create Group Dataset of Percent Change BOLD Map - CP and CnP Trials Across No Report Paradigm

%The purpose of this code is to aggregate the percent change data from
%indivdual subjects.

%Written by: Sharif I. Kronemer
%Date: 5/15/2021

clear

%% Prompts

%Run location
prompt1 = 'Run Location (l/s): ';
run_location = input(prompt1, 's');

%Relevant location
prompt2 = 'NRP Relevant Location (c/q/combine): ';
relevant_location = input(prompt2, 's');

%Center Relevant
if isequal(relevant_location, 'c')
    
    relevant_location = 'Center Relevant';
    relevance_cell = {'Center Relevant'};
    save_loc_name = 'cent';

%Quadrant Relevant
elseif isequal(relevant_location, 'q')
    
    relevant_location = 'Quadrant Relevant';
    relevance_cell = {'Quadrant Relevant'};
    save_loc_name = 'quad';
    
%Combine Center and Quadrant Relevant   
elseif isequal(relevant_location, 'combine')
    
    relevant_location = 'Center and Quadrant Relevant';
	relevance_cell = {'Center Relevant', 'Quadrant Relevant'};
    save_loc_name = 'cent_quad';
    
end

%Stimulus Opacity
prompt3 = 'Stimulus opacity (opaque/thres): ';
stimulus_opacity = input(prompt3, 's');

%% Directories 

%Local
if isequal(run_location, 'l')
    
    %Add directories to path
    addpath(genpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\Extract voxel data\Supplementary functions'));
    addpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\spm12_radiological');
    addpath('S:\Data8\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\Behavioral Analysis')
    
    %Subject direcotry 
    Noreport_subject_folder = 'Y:\HNCT No Report Paradigm\Subject Analysis MRI';
    
%Server
elseif isequal(run_location, 's')
    
    %Add directories to path
    addpath(genpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Extract voxel data/Supplementary functions'));
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/spm12_radiological');
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/Behavioral Analysis')
    addpath('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/fMRI Behavioral Analysis')
    
    %Subject directory
    Noreport_subject_folder = '/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI';
    
    %Save directory
    save_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Data/MRI/BOLD PC';
    
end

%% Group Variables and Subject Lists

%Find all subject folders
Noreport_subject_list = dir(Noreport_subject_folder);
Noreport_subject_list = {Noreport_subject_list.name}';

%Initialized variables
group_CP_BOLD_PC_data = [];
group_CnP_BOLD_PC_data = [];
all_included_subjects_list = {};
CP_total_session_num = 0;
CnP_total_session_num = 0;

%% NO REPORT PARADIGM SUBJECTS - Aggregate percent change over subjects

disp('Aggregate percent change data across subjects')

%Loop over subject
for sub = 1:length(Noreport_subject_list) 
    
    %Select ID
    ID = Noreport_subject_list{sub};
    
    %Initalize subject variables
    subject_CP_BOLD_PC_data = [];
    subject_CnP_BOLD_PC_data = [];
    
    %Loop over relevant locations
    for loc = 1:length(relevance_cell)
        
        %Select relevant location
        location = relevance_cell{loc};
        
        disp(['**Adding No Report Subject ',num2str(ID),' ', location,'**'])
        
        %Rejecting subjects by behavioral performance
        bad_subject_idx = noreport_subject_rejection_by_behavior_MRI(ID, location);
        
        %Check if bad subject by behavior
        if isequal(bad_subject_idx, 1)
            
            disp(['Bad behavior - Skipping ', num2str(ID)])
            continue
            
        end
        
        %Skip subjects where all runs were rejected due to motion
        %rejections (>2mm/1deg) - see motion rejection artifact variable for
        %full list of runs excluded
        
        %Quadrant Relevant
        if ismember(ID,{'632'}) && isequal(location,'Quadrant Relevant')
            
            disp(['All runs rejected by >2mm/1deg motion - Skipping ', num2str(ID)])
            continue
            
        %Center Relevant    
        elseif ismember(ID,{'567','568','648'}) && isequal(location,'Center Relevant')
            
            disp(['All runs rejected by >2mm/1deg motion - Skipping ', num2str(ID)])
            continue
            
        end
        
        %Subject percent change data path 
        if isequal(run_location, 'l')
            
            rootpath = fullfile('Y:\HNCT No Report Paradigm\Subject Analysis MRI',ID,'Perception Task',location,'MRI Session\MRI Analysis\Extracted voxel data\Session percent change\Relevant Stimuli');

        elseif isequal(run_location, 's')
            
            rootpath = fullfile('/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI',ID,'Perception Task',location,'MRI Session/MRI Analysis/Extracted voxel data/Session percent change/Relevant Stimuli');

        end
                    
        %Check if root path exists
        if not(exist(rootpath))

           disp(['Cannot find subject data directory - Skipping ', num2str(ID)])
           continue

        end

        %Load CP percent change data (Note: notstimblinks variable exclude
        %trials with blinks at the time of stimulus onset)
        cd(rootpath)
        load(['extracted_signal_PC_CP_',stimulus_opacity,'_norm_binned_cuts_filtered_nostimblinks.mat'])

        %Add to CP trial count
        CP_total_session_num = CP_total_session_num + size(subject_img_3D_PC_cut,3);
        
        % Create group matrix of time points x voxels x subjects
        if isempty(subject_CP_BOLD_PC_data) 

            subject_CP_BOLD_PC_data(:,:,1)  = subject_img_3D_PC_cut(:,:,:); 

        else

            subject_CP_BOLD_PC_data(:,:,(end+1))  = subject_img_3D_PC_cut(:,:,:);

        end 
     
        %Load CnP percent change data (Note: no stimblinks variable exclude
        %trials with blinks at the time of stimulus onset)
        
        %Threshold opacity trials
        if isequal(stimulus_opacity, 'thres')
            
            cd(rootpath)
            load('extracted_signal_PC_CnP_thres_norm_binned_cuts_filtered_nostimblinks.mat')
                  
            %Add to CnP trial count
            CnP_total_session_num = CnP_total_session_num + size(subject_img_3D_PC_cut,3);
        
        %There are no opaque CnP trials, so replacing this variable with an empty matrix when
        %considering opaque stimuli. 
        else
            
            subject_img_3D_PC_cut = [];
            
        end

        % Create group matrix of time points x voxels x subjects
        if isempty(subject_CnP_BOLD_PC_data) 

            subject_CnP_BOLD_PC_data(:,:,1)  = subject_img_3D_PC_cut(:,:,:); 

        else

            subject_CnP_BOLD_PC_data(:,:,(end+1))  = subject_img_3D_PC_cut(:,:,:);

        end

    end
     
    %% Add subject data to group data matrix
    
    %Check if subject data is empty
    if not(isempty(subject_CP_BOLD_PC_data))
        
        %Add subject ID to subject list
        all_included_subjects_list = [all_included_subjects_list; ID];
    
        %Average within subject if multiple sessions per subject
        if size(subject_CP_BOLD_PC_data,3) > 1
            
            subject_CP_BOLD_PC_data = nanmean(subject_CP_BOLD_PC_data,3);
            subject_CnP_BOLD_PC_data = nanmean(subject_CnP_BOLD_PC_data,3);

        end
        
        % Create group matrix of [time x voxels x subjects]
        if isempty(group_CP_BOLD_PC_data) 

            group_CP_BOLD_PC_data(:,:,1) = subject_CP_BOLD_PC_data; 

        else

            group_CP_BOLD_PC_data(:,:,(end+1)) = subject_CP_BOLD_PC_data;

        end 

        % Create group matrix of [time x voxels x subjects]
        if isempty(group_CnP_BOLD_PC_data) 

            group_CnP_BOLD_PC_data(:,:,1) = subject_CnP_BOLD_PC_data; 

        else

            group_CnP_BOLD_PC_data(:,:,(end+1)) = subject_CnP_BOLD_PC_data;

        end
       
    end
    
end

%% Save Report and No Report Paradigm Dataset Matrices

%Save subject BOLD matrices
cd(save_dir)

%Only saves CP variable for opaque stimuli type
if isequal(stimulus_opacity,'opaque')
    
    save(['Group_PC_BOLD_CP_',save_loc_name,'_',stimulus_opacity,'_data.mat'], 'group_CP_BOLD_PC_data', 'all_included_subjects_list', 'relevant_location', ...
        'CP_total_session_num', 'CnP_total_session_num', '-v7.3')
    
%Saves both CP and CnP files
else
    
    save(['Group_PC_BOLD_CP_CnP_',save_loc_name,'_',stimulus_opacity,'_data.mat'], 'group_CP_BOLD_PC_data', 'group_CnP_BOLD_PC_data', 'all_included_subjects_list', 'relevant_location',...
        'CP_total_session_num', 'CnP_total_session_num', '-v7.3')
    
end
