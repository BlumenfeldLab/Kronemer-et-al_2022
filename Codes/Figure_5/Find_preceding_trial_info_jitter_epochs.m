%% Find Relevant/Opacity Conditions for Preceding Trials in Jitter epochs

%This code is looking at the preceding stimulus for the jitter epochs to
%find if there is a systematic difference between PP and PnP trials. 

%Written by: Sharif I. Kronemer
%Modified: 6/10/2022

clear all

%% Prompts

%Select run location
run_location = 's';

%Run Task Relevant 
relevan_condition = 'rel';

%Stimulus location set
prompt_3 = 'Relevant location [c, q, both]: ';
location_set = input(prompt_3, 's');

%Confidence score threshold
prompt_4 = 'Confidence score threshold [0,0.25,0.5,0.75,1,1.25]: ';
confidence_score_threshold = str2num(input(prompt_4, 's'));

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

%Blink rejection window - stimulus time
%blink_reject_window = [9001:9050]; %[6001:6050];

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
    save_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/MRI Analysis/Percent Change Permutation Stats/Whole Brain_PP_minus_PnP_rel_irrel_cent_quad_jitter_score_thres_0.75_perm_100';
    
    %Movement artifact directory
    motion_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/MRI Analysis/Movement Artifact';
    
end

%% Group Variables and Subject Lists

%Find all subject folders
noreport_subject_list = dir(subject_folder);
noreport_subject_list = {noreport_subject_list.name}';

%Remove non-subject folders
noreport_subject_list([1,2,68,69,70]) = [];

%% Find and Combined PP and PnP Trials Across Subjects
 
disp(['Find PP and PnP trials for score threshold ',num2str(confidence_score_threshold)])

%Number of PP trials
PP_relevant_trials_num = 0;
PnP_relevant_trials_num = 0;
PP_total_num= 0;
PnP_total_num = 0;

PP_CP_num = 0;
PP_CnP_num = 0;

PnP_CP_num = 0;
PnP_CnP_num = 0;

PP_thres_num = 0;
PP_opaque_num = 0;
PP_blank_num = 0;

PnP_thres_num = 0;
PnP_opaque_num = 0;
PnP_blank_num = 0;
        
%Loop over locations sets
for loc = 1:length(location_cell)

    %Select relevant location
    current_location = location_cell{loc}; 
      
    %Load bad runs by excessive motion
    if isequal(current_location, 'Center Relevant') %|| isequal(current_location, 'Quadrant Irrelevant')

        load(fullfile(motion_dir,'MRI_move_artifact_center_rel_2mm_1deg_thresholds.mat'))

    elseif isequal(current_location, 'Quadrant Relevant') %|| isequal(current_location, 'Center Irrelevant')

        load(fullfile(motion_dir,'MRI_move_artifact_quadrant_rel_2mm_1deg_thresholds.mat'))

    end
           
    %Loop over subject
    for sub = 1:length(noreport_subject_list) 

        tic

        %Select current ID
        current_ID = noreport_subject_list{sub};
          
        %% Load Classifier Model Confidence Scores

        %Subject percent change data path 
        if isequal(run_location, 'l')

        elseif isequal(run_location, 's')  

            %Model path
            model_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/Machine Learning Analysis/MRI Session'...
                ,current_location,'Base_boosted_tree_lin_gaus_SVM_stacked_lin_SVM_jitter');
            
        end

        %Load model direcotry 
        cd(model_dir)
        load('Classifier_scores_labels_jitter.mat')
        
        %Store the avg confidence scores into group variable
        group_stacked_SVM_avg_scores = stacked_SVM_avg_scores;

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
        
        %Not processed yet
        if ismember(current_ID,{'568','646'})
           
           disp(['Skipping ',num2str(current_ID)])
           continue
            
        end
        
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
                    
        disp(['**Adding Jitter Trials Subject ',num2str(current_ID),' ', current_location,'**'])      
              
        %% Load Subject Eyelink Data
        
        %Subject data paths
        if isequal(run_location, 'l')
        
        elseif isequal(run_location, 's')

            %Eyelink path
            eyelink_path = fullfile('/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI',current_ID,'Perception Task',current_location,'MRI Session/Pupillometry Analysis');

        end    
              
        %Load eyelink table
        cd(eyelink_path)
        load('eyelinkTable.mat')
        
        %Initialize matrix
        blink_matrix = [];
        relevant_trial = [];
        condition_code = [];
        stim_opacity = [];
        
        %Center relevant or irrelevant
        if isequal(current_location, 'Center Relevant')
                    
            %Loop over eyelink rows
            for row = 1:size(eyelinkTable,1)
                             
                % Relevant first
                if eyelinkTable.RelevantFirst(row) == 1
                    
                    %Blink
                    sub_blink1 = cell2mat(eyelinkTable.BlinkStublinkCenter(row));
                    sub_blink1(1,:) = []; %Remove left eye

                    sub_blink2 = cell2mat(eyelinkTable.BlinkStublinkQuadrant(row));
                    sub_blink2(1,:) = []; %Remove left eye 

                    trial_sub_blink = [sub_blink1; sub_blink2];
                    blink_matrix = [blink_matrix; trial_sub_blink];
                    
                    %Add to relevant trial
                    relevant_trial = [relevant_trial; [1;0]];
                    
                    %Condition code
                    condition_code = [condition_code; [eyelinkTable.ConditionCode(row);NaN]];
                    
                    %Opacity
                    stim_opacity = [stim_opacity; [eyelinkTable.FaceCenterOpacity(row);eyelinkTable.FaceQuadrantOpacity(row)]];
                    
                else
                    
                    %Blink
                    sub_blink1 = cell2mat(eyelinkTable.BlinkStublinkQuadrant(row));
                    sub_blink1(1,:) = []; %Remove left eye

                    sub_blink2 = cell2mat(eyelinkTable.BlinkStublinkCenter(row));
                    sub_blink2(1,:) = []; %Remove left eye 

                    trial_sub_blink = [sub_blink1; sub_blink2];
                    blink_matrix = [blink_matrix; trial_sub_blink];
                    
                    %Add to relevant trial
                    relevant_trial = [relevant_trial; [0;1]];
                    
                    %Condition code
                    condition_code = [condition_code; [NaN;eyelinkTable.ConditionCode(row)]];
                    
                    %Opacity
                    stim_opacity = [stim_opacity; [eyelinkTable.FaceQuadrantOpacity(row);eyelinkTable.FaceCenterOpacity(row)]];
                  
                end

            end
            
        %Quadrant relevant or irrelevant
        elseif isequal(current_location, 'Quadrant Relevant') 

            %Loop over eyelink rows
            for row = 1:size(eyelinkTable,1)
                
                % Relevant first
                if eyelinkTable.RelevantFirst(row) == 1
                    
                    %Blink
                    sub_blink1 = cell2mat(eyelinkTable.BlinkStublinkQuadrant(row));
                    sub_blink1(1,:) = []; %Remove left eye

                    sub_blink2 = cell2mat(eyelinkTable.BlinkStublinkCenter(row));
                    sub_blink2(1,:) = []; %Remove left eye 

                    trial_sub_blink = [sub_blink1; sub_blink2];
                    blink_matrix = [blink_matrix; trial_sub_blink];
                    
                    %Add to relevant trial
                    relevant_trial = [relevant_trial; [1;0]];
                    
                    %Condition code
                    condition_code = [condition_code; [eyelinkTable.ConditionCode(row);NaN]];
                    
                    %Opacity
                    stim_opacity = [stim_opacity; [eyelinkTable.FaceQuadrantOpacity(row);eyelinkTable.FaceCenterOpacity(row)]];
                  
                else
                    
                    %Blink
                    sub_blink1 = cell2mat(eyelinkTable.BlinkStublinkCenter(row));
                    sub_blink1(1,:) = []; %Remove left eye

                    sub_blink2 = cell2mat(eyelinkTable.BlinkStublinkQuadrant(row));
                    sub_blink2(1,:) = []; %Remove left eye 

                    trial_sub_blink = [sub_blink1; sub_blink2];
                    blink_matrix = [blink_matrix; trial_sub_blink];

                    %Add to relevant trial
                    relevant_trial = [relevant_trial; [0;1]];
                                        
                    %Condition code
                    condition_code = [condition_code; [NaN;eyelinkTable.ConditionCode(row)]];
                    
                    %Opacity
                    stim_opacity = [stim_opacity; [eyelinkTable.FaceCenterOpacity(row);eyelinkTable.FaceQuadrantOpacity(row)]];
                                     
                end
                
            end

        end 
                   
        %Extract Run Number
        trial_run_numbers = eyelinkTable.RunNumber;
        
        %% Find the preceding trial info
        
        %Find PP and PnP Trials
        PP_idx = subject_confidence_scores > confidence_score_threshold & stim_opacity > 0 & stim_opacity < 1;
        PnP_idx = subject_confidence_scores < -confidence_score_threshold & stim_opacity > 0 & stim_opacity < 1;
        
        %Count the number of PP and PnP trials
        PP_total_num = PP_total_num + sum(PP_idx);
        PnP_total_num = PnP_total_num + sum(PnP_idx);
        
        %Count the number of relevant trials
        %PP_relevant_trials_num = PP_relevant_trials_num + sum(relevant_trial(PP_idx));
        %PnP_relevant_trials_num = PnP_relevant_trials_num + sum(relevant_trial(PnP_idx));
        
        %Count the condition code
        PP_CP_num = PP_CP_num + sum(condition_code == 1 & PP_idx);
        PP_CnP_num = PP_CnP_num + sum(condition_code == 2 & PP_idx);
        
        PnP_CP_num = PnP_CP_num + sum(condition_code == 1 & PnP_idx);
        PnP_CnP_num = PnP_CnP_num + sum(condition_code == 2 & PnP_idx);
        
        %PP_thres_num = PP_thres_num + sum(stim_opacity > 0 & stim_opacity < 1 & subject_confidence_scores > confidence_score_threshold);
        %PP_opaque_num = PP_opaque_num + sum(stim_opacity > 0.99 & subject_confidence_scores > confidence_score_threshold);
        %PP_blank_num = PP_blank_num + sum(stim_opacity < 0.001 & subject_confidence_scores > confidence_score_threshold);
        
        %PnP_thres_num = PnP_thres_num + sum(stim_opacity > 0 & stim_opacity < 1 & subject_confidence_scores < -confidence_score_threshold);
        %PnP_opaque_num = PnP_opaque_num + sum(stim_opacity > 0.99 & subject_confidence_scores < -confidence_score_threshold);
        %PnP_blank_num = PnP_blank_num + sum(stim_opacity < 0.001 & subject_confidence_scores < -confidence_score_threshold);
                
    end
    
end

%% Calculate Percentage 

PP_opaque_percent = (PP_opaque_num/PP_total_num)*100;
PnP_opaque_percent = (PnP_opaque_num/PnP_total_num)*100;

PP_thres_percent = (PP_thres_num/PP_total_num)*100;
PnP_thres_percent = (PnP_thres_num/PnP_total_num)*100;

PP_blank_percent = (PP_blank_num/PP_total_num)*100;
PnP_blank_percent = (PnP_blank_num/PnP_total_num)*100;

PP_CP_percent = (PP_CP_num/PP_total_num)*100;
PnP_CP_percent = (PnP_CP_num/PnP_total_num)*100;

PP_CnP_percent = (PP_CnP_num/PP_total_num)*100;
PnP_CnP_percent = (PnP_CnP_num/PnP_total_num)*100;

%% Save Data

cd(save_dir)
save trial_type_information.mat PP* PnP*
