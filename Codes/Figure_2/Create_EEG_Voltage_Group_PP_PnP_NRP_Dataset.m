%% Create Group Data of EEG Voltage Data - Report and No Report Paradigm PP and PnP Epochs

%This code will gather all the irrelevant EEG trials across subjects, find the PP and PnP
%based on the selected confidence threshold, average
%within subjects, and then create a group matrix of channel x time x
%subjects. The code also rejects subjects by poor behavioral and rejects
%trials with a blink event during the critical blink window time. 

%Written by: Sharif I. Kronemer
%Date: 5/31/2021
%Modified: 7/12/2021

clear

%% Prompts

%Select run location
prompt_1 = 'Running code local or server [l, s]: ';
run_location = input(prompt_1,'s');

%Run Task Relevant or Irrrelevant Stimuli
% prompt_2 = 'Relevant or irrelevant stimuli [rel, irrel]: ';
% relevant_condition = input(prompt_2, 's');
relevant_condition = 'irrel'; %Currently only functional on irrelevant trials (5/31/2021)

%Stimulus location set
prompt_3 = 'Stimulus location set [c, q, combine]: ';
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

end

%Define opacity value (Note: PP and PnP trials are taken from predictions
%that were made on only threshold opacity trials; opaque trials were not
%considered)
opacity = 'threshold';

%Blink rejection window (-200 +500)
blink_reject_window = [5801:6500]; 

%% Directories 

if isequal(run_location, 's')
    
    %Add group directory to path
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Data/EEG')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/Behavioral Analysis')
    
    %Subject folder directory
    subject_folders = '/mnt/Data28/HNCT No Report Paradigm/Subject Analysis EEG';  
    
    %Save directory
    save_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Data/EEG';
    
elseif isequal(run_location, 'l')
       
end

%% Group Variables and Subject Lists

%Find all subject folders
noreport_subject_list = dir(subject_folders);
noreport_subject_list = {noreport_subject_list.name}';

%Remove non-subject folders
noreport_subject_list([1,2,68,69]) = [];

%Initialized variables
group_PP_EEG_voltage_data = [];
group_PnP_EEG_voltage_data = [];
PP_epochs_subjects_list = {};
PnP_epochs_subjects_list = {};

%Reject trial count
total_rejected_trial_num = 0;

%Total trial count
total_irrelevant_trial_num = 0;
PP_trial_num = 0;
PnP_trial_num = 0;

%% Find and Combined PP and PnP Trials Across Subjects
 
disp(['Finding PP and PnP trials for score threshold ',num2str(confidence_score_threshold)])

tic

%Loop over subject
for sub = 1:length(noreport_subject_list) 

    %Subject SVM irrelevant trial scores
    subject_SVM_scores = [];
    subject_EEG_data = [];
    
    %Select ID
    current_ID = noreport_subject_list{sub};

    %Loop over locations sets
    for loc = 1:length(location_cell)

        %Select relevant location
        current_location = location_cell{loc};
        
        %% Load Subject Trial Index
        
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
            model_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/Machine Learning Analysis/EEG Session'...
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
        bad_subject_idx = noreport_subject_rejection_by_behavior_EEG(current_ID, current_location);
        
        %Check if bad subject by behavior
        if isequal(bad_subject_idx, 1)
            
            disp(['Bad behavior - Skipping ', num2str(current_ID)])
            continue
            
        end
        
        %Exclude subjects without EEG data - Study session not completed
        if ismember(num2str(current_ID), {'623','643','679','680','694'}) && isequal(current_location, 'Center Relevant')

            disp(['Skipping ',num2str(current_ID), ' EEG ', current_location, ' - No Data'])
            continue
            
        elseif ismember(num2str(current_ID), {'623','643','679','680','694'}) && isequal(current_location, 'Quadrant Irrelevant')

            disp(['Skipping ',num2str(current_ID), ' EEG ', current_location, ' - No Data'])
            continue

        elseif ismember(num2str(current_ID), {'581','605','612','653','675'}) && isequal(current_location, 'Quadrant Relevant')

            disp(['Skipping ',num2str(current_ID), ' EEG ', current_location, ' - No Data'])                
            continue
            
        elseif ismember(num2str(current_ID), {'581','605','612','653','675'}) && isequal(current_location, 'Center Irrelevant')

            disp(['Skipping ',num2str(current_ID), ' EEG ', current_location, ' - No Data'])                
            continue
            
        end
        
        disp(['**Adding No Report Subject ',num2str(current_ID),' ', current_location,'**'])

        %% Setup Subject Directories 
        
        %Subject data paths 
        if isequal(run_location, 'l')
        
        elseif isequal(run_location, 's')
            
            %EEG epochs path
            rootpath = fullfile('/mnt/Data28/HNCT No Report Paradigm/Subject Analysis EEG',current_ID,'Perception Task',current_location,'EEG Analysis/Preprocessed Data/Voltage Epochs');

            %Quadrant Relevant and Center Irrelevant    
            if isequal(current_location, 'Quadrant Relevant') || isequal(current_location, 'Center Irrelevant')
                
                %Eyelink path
                eyelink_path = fullfile('/mnt/Data28/HNCT No Report Paradigm/Subject Analysis EEG',current_ID,'Perception Task/Quadrant Relevant/EEG Session/Pupillometry Analysis');
                
                %Bad EEG epochs path
                bad_trial_path = fullfile('/mnt/Data28/HNCT No Report Paradigm/Subject Analysis EEG',current_ID,'Perception Task/Quadrant Relevant/EEG Session/EEG Analysis/Event Times and Identifiers');
            
            %Center Relevant and Quadrant Irrelevant
            elseif isequal(current_location, 'Center Relevant') || isequal(current_location, 'Quadrant Irrelevant')
                
                %Eyelink path
                eyelink_path = fullfile('/mnt/Data28/HNCT No Report Paradigm/Subject Analysis EEG',current_ID,'Perception Task/Center Relevant/EEG Session/Pupillometry Analysis');
                
                %Bad EEG epochs path
                bad_trial_path = fullfile('/mnt/Data28/HNCT No Report Paradigm/Subject Analysis EEG',current_ID,'Perception Task/Center Relevant/EEG Session/EEG Analysis/Event Times and Identifiers');
            
            end
            
        end
                
        %% Load Subject Eyelink Data
        
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
                   
        %Find trials with blink during blink_rejection_window 
        blink_at_stim_trials = find(sum(blink_matrix(:,blink_reject_window),2) > 0);

        %% Load Subject EEG Data        
        
        %Note: 600 Center Relevant has a missing DIN marker in the EEG
        %recordings for the last trial. The Eyelink and irrelevant EEG
        %epcohs already deal with this omission and no adjustments are
        %required here. 

        %Load EEG epoch data - threshold opacity irrelevant/relevant data
        cd(rootpath) 
        
        if isequal(relevant_condition, 'irrel')
        
            load('irrelevant_threshold_4s_epochs.mat')
        
        elseif isequal(relevant_condition, 'rel')
            
           load('relevant_threshold_4s_epochs.mat')
            
        end

        %Remove Run 3 EEG trials due eyelink recording error
        if isequal(current_ID, '652') && isequal(current_location,'Quadrant Irrelevant')

            %Cut trials and indices (Remove run 3: trials 37-54)
            irrelevant_threshold_epochs(:,:,37:54) = [];
            
        elseif isequal(current_ID, '652') && isequal(current_location,'Center Relevant')

            %Cut trials and indices (Remove run 3: trials 37-54)
            relevant_threshold_epochs(:,:,37:54) = [];
            
        end

        %Compare number of trails from EEG and Eyelink data
        if not(isequal(size(irrelevant_threshold_epochs,3), size(eyelinkTable,1)))

            error([current_ID,' Eyelink and EEG do not match - Skipping'])

        end       
        
        %Count total number of irrelevant trials
        total_irrelevant_trial_num = total_irrelevant_trial_num + size(irrelevant_threshold_epochs,3);
        
        %% Rejections by EEG and Eyelink parameters
        
        %Load bad EEG trial index
        cd(bad_trial_path)
        load('Bad_trials_index.mat')
        
        %Update rejection idx - Run 3 EEG bad trials due eyelink recording error
        if isequal(current_ID, '652') && isequal(current_location,'Quadrant Irrelevant')
            
            irrelevant_threshold_bad_epochs_idx([7,8]) = [45,71];
            
        end
        
        %All trials to reject between EEG and Eyelink rejections
        if isequal(relevant_condition,'irrel')
            
            reject_trials = union(blink_at_stim_trials, irrelevant_threshold_bad_epochs_idx);
        
        elseif isequal(relevant_condition,'rel')
            
            reject_trials = union(blink_at_stim_trials, relevant_threshold_bad_epochs_idx);

        end

        %Count the number of bad EEG trials
        total_rejected_trial_num = total_rejected_trial_num + length(reject_trials);
        
        %Reject Trials        
        %Note: Trials are rejected with a blink during critical window or 
        %bad trials according to preprocessing of EEG data
        if isequal(relevant_condition,'irrel')
        
            irrelevant_threshold_epochs(:,:,reject_trials) = []; %EEG data
        
        elseif isequal(relevant_condition,'rel')

            relevant_threshold_epochs(:,:,reject_trials) = []; %EEG data

        end
        
        %SVM scores
        subject_confidence_scores(reject_trials) = [];    
 
        %% Add subject epochs and SVM scores to subject matrix

        % Create subject trials to subject matrix of channels x time x trials
        if isequal(relevant_condition,'irrel')
            
            subject_EEG_data = cat(3, subject_EEG_data, irrelevant_threshold_epochs); 
                
        elseif isequal(relevant_condition,'rel')
            
            subject_EEG_data = cat(3, subject_EEG_data, relevant_threshold_epochs); 
        
        end

        %Cut subject SVM scores from group matrix and add to subject score
        %matrix
        subject_SVM_scores = cat(1,subject_SVM_scores,subject_confidence_scores);
        
    end
    
    %% Confirm equal trial number between EEG data and SVM scores
    
    %Check if any EEG data
    if isempty(subject_EEG_data)
    
       disp(['No EEG data for subject ',num2str(current_ID),' - Skipping from group addition'])
       continue
    
    elseif not(isequal(length(subject_SVM_scores), size(subject_EEG_data,3)))
        
       error(['Trial mismatch for subject ',num2str(current_ID)]) 
        
    end
    
    %% Find the Predicted Trials For Selected Confidence Score Threshold

    %Extract the predicted class trials
    PP_subject_EEG = subject_EEG_data(:,:,subject_SVM_scores > confidence_score_threshold);
    PnP_subject_EEG = subject_EEG_data(:,:,subject_SVM_scores < -confidence_score_threshold);
    
    %Count the number of PP and PnP
    PP_trial_num = PP_trial_num + size(PP_subject_EEG,3);
    PnP_trial_num = PnP_trial_num + size(PnP_subject_EEG,3);

    %% Store Average Predicted Trials in Group Matrix

    %Add subject PP data to group 
    if ~isempty(PP_subject_EEG)   

        %Add subject ID to subject list
        PP_epochs_subjects_list = [PP_epochs_subjects_list; current_ID];

        %Average within subject trials if multiple trials per subject
        if size(PP_subject_EEG,3) > 1

            PP_subject_EEG = nanmean(PP_subject_EEG,3);

        end
        
    else
       
       %If no data create an NAN matrix 
       PP_subject_EEG = nan(257,4001);
        
    end
    
    % Create group matrix of [channel x time x subjects]
    if isempty(group_PP_EEG_voltage_data) 

        group_PP_EEG_voltage_data(:,:,1) = PP_subject_EEG; 

    else

        group_PP_EEG_voltage_data(:,:,(end+1)) = PP_subject_EEG;

    end 

    %Add subject PnP data to group 
    if ~isempty(PnP_subject_EEG)   

        %Add subject ID to subject list
        PnP_epochs_subjects_list = [PnP_epochs_subjects_list; current_ID];

        %Average within subject trials if multiple trials per subject
        if size(PnP_subject_EEG,3) > 1

            PnP_subject_EEG = nanmean(PnP_subject_EEG,3);

        end

    else
        
       %If no data create an NAN matrix 
       PnP_subject_EEG = nan(257,4001);
           
    end
    
    % Create group matrix of [channel x time x subjects]
    if isempty(group_PnP_EEG_voltage_data) 

        group_PnP_EEG_voltage_data(:,:,1) = PnP_subject_EEG; 

    else

        group_PnP_EEG_voltage_data(:,:,(end+1)) = PnP_subject_EEG;

    end              
      
end

toc

%% Save Group Data

%Data type folder
save_folder_name = 'Voltage';

%Save group EEG data
cd(fullfile(save_dir,save_folder_name))

save(['Group_EEG_voltage_',condition_name,'_PP_PnP_',rel_save_name,'_',opacity,'_score_thres_',num2str(confidence_score_threshold),'_data.mat'],...
    'PP_epochs_subjects_list','PnP_epochs_subjects_list','group_PP_EEG_voltage_data','group_PnP_EEG_voltage_data',...
    'PP_trial_num','PnP_trial_num','total_rejected_trial_num','total_irrelevant_trial_num','-v7.3')    
