%% Find PP and PnP Trials Across No Report Paradigm

%The purpose of this code is to aggregate the percent change data from
%indivdual subjects.

%Written by: Sharif I. Kronemer
%Date: 6/1/2021
%Modified: 9/24/2021

clear

%% Prompts

%Select run location
% prompt_1 = 'Running code local or server [l, s]: ';
% run_location = input(prompt_1,'s');
run_location = 's';

%Run Task Relevant or Irrrelevant Stimuli
% prompt_2 = 'Relevant or irrelevant stimuli [rel, irrel]: ';
% relevant_condition = input(prompt_2, 's');
relevant_condition = 'irrel'; %Currently only functional on irrelevant trials (5/31/2021)

%Stimulus location set
prompt_3 = 'Stimulus irrelevant location set [c, q, both]: ';
location_set = input(prompt_3, 's');

%Confidence score threshold
prompt_4 = 'Confidence score threshold [0,0.25,0.5,0.75,1,1.25]: ';
confidence_score_threshold = str2num(input(prompt_4, 's'));

%Relevant stimuli
if isequal(relevant_condition, 'rel')
    
    %Condition name
    condition_name = 'relevant';
    
    %Center Relevant
    if isequal(location_set, 'c')

        location_cell = {'Center Relevant'};

    %Quadrant Relevant
    elseif isequal(location_set, 'q')

        location_cell = {'Quadrant Relevant'};

    %Both Center and Quadrant Relevant   
    elseif isequal(location_set, 'both')

        location_cell = {'Center Relevant','Quadrant Relevant'};

    end

%Irrelevant stimuli
elseif isequal(relevant_condition, 'irrel')
    
    %Condition name
    condition_name = 'irrelevant';
    
    %Center Irrelevant
    if isequal(location_set, 'c')

        location_cell = {'Center Irrelevant'};

    %Quadrant Irrelevant
    elseif isequal(location_set, 'q')

        location_cell = {'Quadrant Irrelevant'};

    %Both Center and Quadrant Irelevant   
    elseif isequal(location_set, 'both')

        location_cell = {'Center Irrelevant','Quadrant Irrelevant'};

    end

end

%Define opacity value (Note: PP and PnP trials are taken from predictions
%that were made on only threshold opacity trials; opaque trials were not
%considered)
opacity = 'threshold';

%Blink rejection window - stimulus time
blink_reject_window = [6001:6050]; 

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
    subject_folder = '/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI';
    
    %Save directory
    save_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Data/MRI';
    
    %Movement artifact directory
    motion_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/MRI Analysis/Movement Artifact';
    
end

%% Group Variables and Subject Lists

%Find all subject folders
noreport_subject_list = dir(subject_folder);
noreport_subject_list = {noreport_subject_list.name}';

%Remove non-subject folders
noreport_subject_list([1,2,68,69,70]) = [];

%Total trial count
PP_trial_num = 0;
PnP_trial_num = 0;

%% Find and Combined PP and PnP Trials Across Subjects
 
disp(['Find PP and PnP trials for score threshold ',num2str(confidence_score_threshold)])

%Loop over locations sets
for loc = 1:length(location_cell)

    %Select relevant location
    current_location = location_cell{loc}; 
      
    %Load bad runs by excessive motion
    if isequal(current_location, 'Center Relevant') || isequal(current_location, 'Quadrant Irrelevant')

        load(fullfile(motion_dir,'MRI_move_artifact_center_rel_2mm_1deg_thresholds.mat'))

    elseif isequal(current_location, 'Quadrant Relevant') || isequal(current_location, 'Center Irrelevant')

        load(fullfile(motion_dir,'MRI_move_artifact_quadrant_rel_2mm_1deg_thresholds.mat'))

    end
           
    %Loop over subject
    for sub = 1:length(noreport_subject_list) 

        tic
        
        %Initialize Subject Variables
        PP_BOLD = [];
        PnP_BOLD = [];
        PP_BOLD_time = [];
        PnP_BOLD_time = [];

        %Select current ID
        current_ID = noreport_subject_list{sub};
    
        %% Load Subject Trial Index
        
        %Note: The subject index is a list of subject IDs for each trial
        %among the subject aggregated trials
        
        %Location of trial index
        cd(save_dir)
        
        %Relevant condition
        if isequal(relevant_condition, 'rel')
        
            %Quadrant Relevant
            if isequal(current_location, 'Quadrant Relevant')
                
                load('subject_trial_index_quadrant_relevant.mat')
            
            %Center Relevant
            elseif isequal(current_location, 'Center Relevant')
            
                load('subject_trial_index_center_relevant.mat')
            
            end
        
        %Irrelevant condition    
        elseif isequal(relevant_condition, 'irrel')
            
            %Quadrant Irrelevant
            if isequal(current_location, 'Quadrant Irrelevant')
                
                load('subject_trial_index_quadrant_irrelevant.mat')
            
            %Center Irrelevant
            elseif isequal(current_location, 'Center Irrelevant')
            
                load('subject_trial_index_center_irrelevant.mat')
            
            end
            
        end
    
         %% Load Classifier Model Confidence Scores

        %Subject percent change data path 
        if isequal(run_location, 'l')

        elseif isequal(run_location, 's')  

            %Model path
            model_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/Machine Learning Analysis/MRI Session'...
                ,current_location,'Base_boosted_tree_lin_gaus_SVM_stacked_lin_SVM');
            
        end

        %Load model direcotry 
        cd(model_dir)

        %Relevant data model
        if isequal(relevant_condition, 'rel')

            load('Classifier_models_scores_labels.mat')

            %Store the avg confidence scores into group variable
            group_stacked_SVM_avg_scores = stacked_scores;

        %Irrelevant data model
        elseif isequal(relevant_condition, 'irrel')

            load('Classifier_noreport_scores_labels.mat')

            %Store the avg confidence scores into group variable
            group_stacked_SVM_avg_scores = stacked_SVM_avg_scores;

        end
        
        %Extract subject specific confidence scores using subject_trial_index
        subject_confidence_scores = group_stacked_SVM_avg_scores(find(strcmp(subject_trial_index,current_ID)));
     
        %% Subject Exclusions
        
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
              
        %% Load Subject Eyelink Data
        
        %Subject data paths
        if isequal(run_location, 'l')
        
        elseif isequal(run_location, 's')
            
            %Quadrant Relevant and Center Irrelevant    
            if isequal(current_location, 'Quadrant Relevant') || isequal(current_location, 'Center Irrelevant')
                
                %Eyelink path
                eyelink_path = fullfile('/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI',current_ID,'Perception Task/Quadrant Relevant/MRI Session/Pupillometry Analysis');

                %Center Relevant and Quadrant Irrelevant
            elseif isequal(current_location, 'Center Relevant') || isequal(current_location, 'Quadrant Irrelevant')
                
                %Eyelink path
                eyelink_path = fullfile('/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI',current_ID,'Perception Task/Center Relevant/MRI Session/Pupillometry Analysis');
                
            end
            
        end    
              
        %Load eyelink table
        cd(eyelink_path)
        load('eyelinkTable.mat')

        %Center relevant or irrelevant
        if isequal(current_location, 'Center Relevant') || isequal(current_location, 'Center Irrelevant')

            %Keep only threshold trials
            eyelinkTable(eyelinkTable.FaceCenterOpacity == 0,:) = []; %Blank trials
            eyelinkTable(eyelinkTable.FaceCenterOpacity == 1,:) = []; %Opaque trials
            eyelinkTable(eyelinkTable.FaceCenterOpacity > 1,:) = []; %For 618 & 631 MRI that has a run of opacity greater that 1

            %Create blink matrix
            
            %Convert to matrix and remove left eye
            blink_matrix = cell2mat(eyelinkTable.BlinkStublinkCenter);
            blink_matrix(1:2:size(eyelinkTable,1)*2,:) = [];   
            
        %Quadrant relevant or irrelevant
        elseif isequal(current_location, 'Quadrant Relevant') || isequal(current_location, 'Quadrant Irrelevant')

            %Keep only threshold trials
            eyelinkTable(eyelinkTable.FaceQuadrantOpacity == 0,:) = []; %Blank trials
            eyelinkTable(eyelinkTable.FaceQuadrantOpacity == 1,:) = []; %Opaque trials
            eyelinkTable(eyelinkTable.FaceQuadrantOpacity > 1,:) = []; %For 618 & 631 MRI that has a run of opacity greater that 1
            
            %Create blink matrix
            
            %Convert to matrix and remove left eye
            blink_matrix = cell2mat(eyelinkTable.BlinkStublinkQuadrant);
            blink_matrix(1:2:size(eyelinkTable,1)*2,:) = [];   

        end 
                   
        %Extract Run Number
        trial_run_numbers = eyelinkTable.RunNumber;
        
        %% Load Irrelevant Face Event Times
    
        %Define event directory
        if isequal(run_location,'l')
            
        elseif isequal(run_location,'s')
        
            event_dir = ['/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI/',current_ID,'/Perception Task/',current_location,'/MRI Analysis/Event Times'];       
        
        end
        
        %Load event times
        cd(event_dir)
        load MRI_irrelevant_face_event_times.mat threshold*
    
        %% Load MRI Data per Run
        
        %Find subject runs
        sub_runs = unique(trial_run_numbers);
    
        %Loop over runs
        for run = 1:length(sub_runs)

            %Find current run from sub_run index
            current_run = sub_runs(run);

            %Check if Run is excluded by excessive movement
            if ismember([current_ID,' Run_',num2str(current_run)],num_bad_runs_by_mm_deg)
                
                disp(['Skipping Run ',num2str(current_run),' - Excessive movement'])
                continue
            
            else
                
                disp(['Running Run ',num2str(current_run)])

            end
        
            %Define Session and Run BOLD Directories
            if isequal(run_location,'l')

            elseif isequal(run_location,'s')

                %Irrelevant stimulus condition
                if isequal(relevant_condition, 'irrel')
                
                    %Run directory
                    BOLD_PC_run_data = fullfile('/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI',current_ID,'Perception Task',current_location,...
                        'MRI Analysis/Extracted voxel data',['Run_',num2str(current_run)],'axialSlices/Irrelevant Stimuli');

                    %Session directory 
                    BOLD_PC_session_data = fullfile('/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI',current_ID,'Perception Task',current_location,...
                        'MRI Analysis/Extracted voxel data/Session percent change/Irrelevant Stimuli');

                end
                
            end            
            
            %Load irrelevant threshold stimuli BOLD
            load(fullfile(BOLD_PC_run_data,'extracted_signal_PC_threshold_norm_binned_cuts_filtered.mat'))

            %Cut SVM scores and blink matrix to current run trials
            run_confidence_scores = subject_confidence_scores(trial_run_numbers == current_run);
            run_blink_matrix = blink_matrix(trial_run_numbers == current_run,:);

            %Find trials with blinks during blink_rejection_window for current run number
            blink_at_stim_trials = find(sum(run_blink_matrix(:,blink_reject_window),2) > 0);                       
            
            %Remove trials with blinks during rejection window from BOLD,
            %times, and scores matrices
            run_confidence_scores(blink_at_stim_trials) = [];
            eval(['threshold_face_times_run_',num2str(current_run),'(blink_at_stim_trials) = [];'])
            times_img_3D_PC_cut(:,blink_at_stim_trials) = [];
            img_3D_PC_cut(:,:,blink_at_stim_trials) = [];
            
            %Confirm appropriate index and data trial numbers
            if not(isequal(size(run_confidence_scores,1), size(img_3D_PC_cut,3), size(times_img_3D_PC_cut,2),...
                    size(eval(['threshold_face_times_run_',num2str(current_run)]),1)))

                error('**Number of trials mismatch with BOLD and index**')

            end
            
            %Find predicted perceived and not perceived trials for
            %specified SVM score threshold
            PP_BOLD = cat(3, PP_BOLD, img_3D_PC_cut(:,:,run_confidence_scores > confidence_score_threshold));
            PnP_BOLD = cat(3, PnP_BOLD, img_3D_PC_cut(:,:,run_confidence_scores < -confidence_score_threshold));
            PP_BOLD_time = cat(2, PP_BOLD_time, times_img_3D_PC_cut(:,run_confidence_scores > confidence_score_threshold));
            PnP_BOLD_time = cat(2, PnP_BOLD_time, times_img_3D_PC_cut(:,run_confidence_scores < -confidence_score_threshold));

            %Clear BOLD data variable
            clearvars img_3D_PC_cut
            
            %Find PP and PnP trials event times for each run
            eval(['PP_event_times_run_',num2str(current_run),' = threshold_face_times_run_',num2str(current_run),'(run_confidence_scores > confidence_score_threshold);'])
            eval(['PnP_event_times_run_',num2str(current_run),' = threshold_face_times_run_',num2str(current_run),'(run_confidence_scores < -confidence_score_threshold);'])

        end    
        
        toc
        
        %Save the PP PnP event times
        cd(event_dir) 
        save(['PP_PnP_irrelevant_face_event_times_confidence_score_',num2str(confidence_score_threshold),'.mat'],'PP_event*','PnP_event*')
    
        %% Bin and recut BOLD epochs data

        %Generate vector for each bin each vector size is nvoxels*1.
        pre_face_onset = 20; %Seconds pre task event onset
        post_face_onset = 20; %Seconds post task event onset

        %Define the different data types to bin
        data_types = {'PP_BOLD','PnP_BOLD'};

        %Loop over types of epochs
        for type = 1:length(data_types)

            %Select trial type
            current_data_type = data_types{type}; 

           %Check of PP matrix is empty
            if not(isempty(eval(current_data_type))) && strcmp(current_data_type,'PP_BOLD')
                
                %Number of PP trials
                PP_trial_num = size(PP_BOLD,3);

                disp(['Binning ',current_data_type, ' in progress'])
            
            %Check if PnP matrix is empty
            elseif not(isempty(eval(current_data_type))) && strcmp(current_data_type,'PnP_BOLD')

                %Number of PnP trials
                PnP_trial_num = size(PnP_BOLD,3);
                               
                disp(['Binning ',current_data_type, ' in progress'])

            else

               disp(['No MRI data for subject ',num2str(current_ID),' - Skipping from binning and create empty matrix'])
               
               %Create empty variable
               eval([current_data_type,'_subject_img_3D_PC_cut = [];'])
               continue       

            end    

            %Initialize bins
            for bin = (-pre_face_onset):1:post_face_onset

                eval(['bin_time_' num2str(bin+pre_face_onset) ' = [];'])

            end

            %Rename time variable
            current_times = eval([current_data_type,'_time']);

            %Find the number of faces and loop over faces
            for face_index = 1:size(current_times,2)

                %Find the number to timepoints/volumes and loop over time
                for vol_time_index = 1:size(current_times,1) 

                    %Find the specific volume time for a particular trial
                    vol_time = current_times(vol_time_index,face_index);

                    %Loop over bins/time
                    for bin = (-pre_face_onset):1:post_face_onset

                        %Find the volume times nearest bin
                        if (vol_time>=bin) && (vol_time<(bin+1)) %new bining method 

                            %Store XYZ MR data into specific bin matrix
                            eval(['bin_time_' num2str(bin+pre_face_onset) '(:,end+1) = ',current_data_type,'(vol_time_index,:,face_index);'])

                        end

                    end

                end

            end

            %Clear BOLD PP/PnP variable for memory space
            eval(['clearvars ',  current_data_type]) 

            %Get average for each bin. 
            subject_img_3D_PC_cut = []; %dimensions are #voxels * time

            %Average XYZ signal over faces
            for bin = (-pre_face_onset):1:post_face_onset

                %Rare cases when bin_time_X is empty, fill with nan before
                %adding to the session level matrix subject_img_3D_PC_cut
                if isempty(eval(['bin_time_',num2str(bin+pre_face_onset)]))

                    disp(['Replacing empty bin_time_', num2str(bin+pre_face_onset),' with nan'])

                    %Enter a voxel matrix of NaNs to replace empty bin time
                    eval(['subject_img_3D_PC_cut(end+1,:) = nan(1,902629);'])

                else

                    %Enter average voxel values over faces 
                    %Dimensions of subject_img_3D_PC_cut are time*#voxels 
                    eval(['subject_img_3D_PC_cut(end+1,:) = nanmean(bin_time_' num2str(bin+pre_face_onset) ',2);'])

                end

            end

            %Rename binned BOLD variable
            eval([current_data_type,'_subject_img_3D_PC_cut = subject_img_3D_PC_cut;'])

            %Clear bins for memory space
            clearvars bin*

        end

        %Save subject predicted PP/PnP binned BOLD data and trial number
        cd(BOLD_PC_session_data)
        save(['extracted_signal_PC_predicted_threshold_opacity_score_',num2str(confidence_score_threshold),'.mat'], ...
            'PP_BOLD_subject_img_3D_PC_cut', 'PnP_BOLD_subject_img_3D_PC_cut', 'PP_trial_num', 'PnP_trial_num', '-v7.3')
 
        toc
        
        %Clear variables before continuing
        clearvars PP_BOLD* PnP_BOLD*
        
    end  
    
end
