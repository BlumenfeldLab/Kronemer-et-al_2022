%% Create Group Dataset of Percent Change BOLD Map - PP and PnP Trials Across No Report Paradigm Jitter

%The purpose of this code is to aggregate the percent change data from
%indivdual subjects.

%Written by: Sharif I. Kronemer
%Date: 6/2/2021
%Last modified: 6/27/2022

clear

%% Prompts

%Run location
prompt1 = 'Run Location (l/s): ';
run_location = input(prompt1, 's');

%Run Task Relevant or Irrrelevant Stimuli
relevant_condition = 'rel'; 

%Relevant location
prompt3 = 'NRP Relevant Location (c/q/combine): ';
location_set = input(prompt3, 's');

%Confidence score threshold
prompt_4 = 'Confidence score threshold [0,0.25,0.5,0.75,1,1.25]: ';
confidence_score_threshold = str2num(input(prompt_4, 's'));
    
%Condition name
condition_name = 'rel_irrel';

%Center Relevant
if isequal(location_set, 'c')

    location_cell = {'Center Relevant'};
    save_loc_name = 'cent';

%Quadrant Relevant
elseif isequal(location_set, 'q')

    location_cell = {'Quadrant Relevant'};
    save_loc_name = 'quad';

%Both Center and Quadrant Relevant   
elseif isequal(location_set, 'combine')

    location_cell = {'Center Relevant','Quadrant Relevant'};
    save_loc_name = 'cent_quad';

end

%% Directories 

%Local
if isequal(run_location, 'l')
    
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
noreport_subject_list = dir(Noreport_subject_folder);
noreport_subject_list = {noreport_subject_list.name}';

%Remove non-subject folders
noreport_subject_list([1,2,68,69,70]) = [];

%Initialized variables
group_PP_BOLD_PC_data = [];
group_PnP_BOLD_PC_data = [];
PP_epochs_subjects_list = {};
PnP_epochs_subjects_list = {};

%Start trial count variables
total_PP_trial_num = 0;
total_PnP_trial_num = 0;
total_PP_unfilled_trial_num = 0;
total_PnP_unfilled_trial_num = 0;

%% NO REPORT PARADIGM SUBJECTS - Aggregate percent change over subjects

tic

disp('Aggregate percent change data across subjects')

%Loop over subject
for sub = 1:length(noreport_subject_list) 
    
    %Select ID
    current_ID = noreport_subject_list{sub};
    
    %Initalize subject variables
    subject_PP_BOLD_PC_data = [];
    subject_PnP_BOLD_PC_data = [];
    
    %Loop over relevant locations
    for loc = 1:length(location_cell)
        
        %Select relevant location
        current_location = location_cell{loc};
               
        %Rejecting subjects by behavioral performance
        bad_subject_idx = noreport_subject_rejection_by_behavior_MRI(current_ID, current_location);
        
        %Check if bad subject by behavior
        if isequal(bad_subject_idx, 1)
            
            disp(['Bad behavior - Skipping ', num2str(current_ID)])
            continue
            
        end
        
        %Skip subjects where all runs were rejected due to motion
        %rejections (>2mm/1deg) - see motion rejection artifact variable for
        %full list of runs excluded
        
        %Quadrant Relevant
        if ismember(current_ID,{'632'}) && ismember(current_location,{'Quadrant Relevant','Center Irrelevant'})
            
            disp(['All runs rejected by >2mm/1deg motion - Skipping ', num2str(current_ID)])
            continue
            
        %Center Relevant    
        elseif ismember(current_ID,{'567','568','648'}) && ismember(current_location,{'Center Relevant','Quadrant Irrelevant'})
            
            disp(['All runs rejected by >2mm/1deg motion - Skipping ', num2str(current_ID)])
            continue
            
        end
        
        %Exclude subjects without MRI data - Study session not completed
        if ismember(num2str(current_ID), {'579','610','623','643','679'}) && ismember(current_location, {'Center Relevant','Quadrant Irrelevant'})

            disp(['Skipping ',num2str(current_ID), ' MRI ', current_location, ' - No Data'])
            continue
            
        elseif ismember(num2str(current_ID), {'581','587','605','612','653'}) && ismember(current_location, {'Quadrant Relevant','Center Irrelevant'})

            disp(['Skipping ',num2str(current_ID), ' MRI ', current_location, ' - No Data'])                
            continue
            
        end
                    
        disp(['**Adding No Report Subject ',num2str(current_ID),' ', current_location,'**'])      
        
        %% Load Subject PP and PnP BOLD 
        
        %BOLD directory
        BOLD_dir = fullfile('/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI',current_ID,'Perception Task',current_location,...
            'MRI Session/MRI Analysis/Extracted voxel data/Session percent change/Relevant Stimuli');
            
        %Load PP/PnP BOLD
        cd(BOLD_dir)
        
        try
            
            %load(['extracted_signal_PC_predicted_jitter_score_',num2str(confidence_score_threshold),'.mat'])
            load(['extracted_signal_PC_predicted_jitter_matched_score_',num2str(confidence_score_threshold),'.mat'])

        catch
        
            disp(['Skipping subject ',num2str(current_ID)])
            continue
        
        end
        
        %Check if PP BOLD data is empty before adding it to group
        if isempty(PP_BOLD_subject_img_3D_PC_cut) 
            
            disp(['No PP BOLD data - not adding to group'])
        
        else
            
            %Count the number of PP trials
            total_PP_trial_num = total_PP_trial_num + PP_trial_num;
            total_PP_unfilled_trial_num = total_PP_unfilled_trial_num + PP_unfilled_num;

            % Create group matrix of time points x voxels x subjects
            if isempty(subject_PP_BOLD_PC_data) 

                subject_PP_BOLD_PC_data(:,:,1)  = PP_BOLD_subject_img_3D_PC_cut(:,:,:); 

            else

                subject_PP_BOLD_PC_data(:,:,(end+1)) = PP_BOLD_subject_img_3D_PC_cut(:,:,:);

            end 
     
        end
        
        %Check if PnP BOLD data is empty before addign it to group
        if isempty(PnP_BOLD_subject_img_3D_PC_cut) 
            
            disp(['No PnP BOLD data - not adding to group'])
        
        else
            
            %Count the number of PP trials
            total_PnP_trial_num = total_PnP_trial_num + PnP_trial_num;
            total_PnP_unfilled_trial_num = total_PnP_unfilled_trial_num + PnP_unfilled_num;
            
            % Create group matrix of time points x voxels x subjects
            if isempty(subject_PnP_BOLD_PC_data) 

                subject_PnP_BOLD_PC_data(:,:,1)  = PnP_BOLD_subject_img_3D_PC_cut(:,:,:); 

            else

                subject_PnP_BOLD_PC_data(:,:,(end+1)) = PnP_BOLD_subject_img_3D_PC_cut(:,:,:);

            end 
     
        end
        
        %Clear current BOLD data
        clearvars PP_BOLD* PnP_BOLD* PP_trial_num PnP_trial_num

    end
     
    %% Add subject data to group data matrix
    
    %Check if subject PP data is empty
    if not(isempty(subject_PP_BOLD_PC_data))
        
        %Add subject ID to PP subject list
        PP_epochs_subjects_list = [PP_epochs_subjects_list; {current_ID}];
    
        %Average within subject if multiple sessions per subject
        if size(subject_PP_BOLD_PC_data,3) > 1
            
            subject_PP_BOLD_PC_data = nanmean(subject_PP_BOLD_PC_data,3);

        end
        
        % Create group matrix of [time x voxels x subjects]
        if isempty(group_PP_BOLD_PC_data) 

            group_PP_BOLD_PC_data(:,:,1) = subject_PP_BOLD_PC_data; 

        else

            group_PP_BOLD_PC_data(:,:,(end+1)) = subject_PP_BOLD_PC_data;

        end
        
    end
       
    %Check if subject PnP data is empty
    if not(isempty(subject_PnP_BOLD_PC_data))
        
        %Add subject ID to PP subject list
        PnP_epochs_subjects_list = [PnP_epochs_subjects_list; {current_ID}];
    
        %Average within subject if multiple sessions per subject
        if size(subject_PnP_BOLD_PC_data,3) > 1
            
            subject_PnP_BOLD_PC_data = nanmean(subject_PnP_BOLD_PC_data,3);

        end
        
        % Create group matrix of [time x voxels x subjects]
        if isempty(group_PnP_BOLD_PC_data) 

            group_PnP_BOLD_PC_data(:,:,1) = subject_PnP_BOLD_PC_data; 

        else

            group_PnP_BOLD_PC_data(:,:,(end+1)) = subject_PnP_BOLD_PC_data;

        end
        
    end
    
end

toc

%% Save Report and No Report Paradigm Dataset Matrices

%Save subject BOLD matrices
cd(save_dir)
%save(['Group_PC_BOLD_',condition_name,'_PP_PnP_',save_loc_name,'_jitter_score_thres_',num2str(confidence_score_threshold),'_data.mat'],...
%    'group_PP_BOLD_PC_data', 'group_PnP_BOLD_PC_data', 'PP_epochs_subjects_list', 'PnP_epochs_subjects_list', 'total_PP_trial_num', 'total_PnP_trial_num', '-v7.3')
save(['Group_PC_BOLD_',condition_name,'_PP_PnP_',save_loc_name,'_jitter_matched_score_thres_',num2str(confidence_score_threshold),'_data.mat'],...
    'group_PP_BOLD_PC_data', 'group_PnP_BOLD_PC_data', 'PP_epochs_subjects_list', 'PnP_epochs_subjects_list', 'total_PP_trial_num', 'total_PnP_trial_num',...
    'total_PP_unfilled_trial_num', 'total_PnP_unfilled_trial_num','-v7.3')
                                