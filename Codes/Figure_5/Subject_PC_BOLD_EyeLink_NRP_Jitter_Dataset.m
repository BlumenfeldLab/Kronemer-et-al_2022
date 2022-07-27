 %% Find EyeLink Trials Across No Report Paradigm Jitter Trials - Pupil, Blink, and Microsaccade

%This code will look for the pupil, blink, and microsaccade selected trials
%from the jitter trials and create subject BOLD epochs for each. 

%Written by: Sharif I. Kronemer
%Modified: 6/28/2022

clear

%% Prompts

%Select run location
run_location = 's';

%Stimulus location set
prompt_3 = 'Relevant location [c, q, both]: ';
location_set = input(prompt_3, 's');

%Data type
prompt_5 = 'Data type [pupil, blink, microsac]: ';
eyelink_type = input(prompt_5,'s');

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

%Blink rejection window - stimulus time - Original is 6001:6050 for the stimulus time
blink_reject_window = [9001:9050];

%Match sample size to main figure dataset (yes = 1; no = 0)
%main_fig_PP_num = 1785;
%main_fig_PnP_num = 1504;
match_sample_size = 1;

%PP/PnP matching ratio
if isequal(eyelink_type,'pupil')

    %all_jitter_PP_num = 5106;
    %all_jitter_PnP_num = 9144;
    PP_ratio = 0.4;%0.28;
    PnP_ratio = 0.11;%0.08;
    
elseif isequal(eyelink_type,'blink')

    %all_jitter_PP_num = 3856;
    %all_jitter_PnP_num = 8913;
    PP_ratio = 0.4;%0.24;
    PnP_ratio = 0.12;%0.09;

elseif isequal(eyelink_type,'microsac')

    %all_jitter_PP_num = 8085;
    %all_jitter_PnP_num = 7721;
    PP_ratio = 0.15;
    PnP_ratio = 0.25;

end

%PP_ratio = 0.11 %main_fig_PP_num/all_jitter_PP_num;
%PnP_ratio = 0.18 %main_fig_PnP_num/all_jitter_PnP_num;

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
    
    %Movement artifact directory
    motion_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/MRI Analysis/Movement Artifact';
    
    %Trial index dir
    index_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/EyeLink Analysis/EyeLink Event Times/MRI/Jitter Epochs';

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
all_keep_PP_num = 0;
all_keep_PnP_num = 0;

%% Find and Combined PP and PnP Trials Across Subjects
 
%Loop over locations sets
for loc = 1:length(location_cell)

    %Select relevant location
    current_location = location_cell{loc}; 
      
    %Load bad runs by excessive motion
    if isequal(current_location, 'Center Relevant') 

        load(fullfile(motion_dir,'MRI_move_artifact_center_rel_2mm_1deg_thresholds.mat'))

    elseif isequal(current_location, 'Quadrant Relevant') 

        load(fullfile(motion_dir,'MRI_move_artifact_quadrant_rel_2mm_1deg_thresholds.mat'))

    end
           
    %Load jitter trial index/list
    cd(index_dir)

    if isequal(current_location, 'Center Relevant')

        load('eyelink_center_rel_selected_jitter_trials.mat')

    elseif isequal(current_location, 'Quadrant Relevant')

        load('eyelink_quadrant_rel_selected_jitter_trials.mat')

    end

    % Select the desired trial list
    if strcmp(eyelink_type, 'pupil')

        trial_IDX = pupil_trial_condition_idx;

    elseif strcmp(eyelink_type, 'blink')

        trial_IDX = blink_trial_condition_idx;

    elseif strcmp(eyelink_type, 'microsac')

        trial_IDX = microsaccade_trial_condition_idx;

    end

    %Store trial_IDX to group variable
    group_trial_IDX = trial_IDX;
    
    %Loop over subject
    for sub = 1:length(noreport_subject_list) %[32,61]%16

        tic
        
        %Initialize Subject Variables
        PP_BOLD = [];
        PnP_BOLD = [];
        PP_BOLD_time = [];
        PnP_BOLD_time = [];

        %Select current ID
        current_ID = noreport_subject_list{sub};
          
        %% Load EyeLink selected Trials

        %Extract subject specific confidence scores using subject_trial_index
        subject_eyelink_idx = group_trial_IDX(find(strcmp(subject_trial_index,current_ID)));

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
        current_subject_opacity = [];
        
        %Center relevant or irrelevant
        if isequal(current_location, 'Center Relevant') 
            
            %Loop over eyelink rows
            for row = 1:size(eyelinkTable,1)

                % Relevant first
                if eyelinkTable.RelevantFirst(row) == 1

                    sub_blink1 = cell2mat(eyelinkTable.BlinkStublinkCenter(row));
                    sub_blink1(1,:) = []; %Remove left eye

                    sub_blink2 = cell2mat(eyelinkTable.BlinkStublinkQuadrant(row));
                    sub_blink2(1,:) = []; %Remove left eye 
                    
                    sub_opacity1 = eyelinkTable.FaceCenterOpacity(row);
                    sub_opacity2 = eyelinkTable.FaceQuadrantOpacity(row);

                    trial_sub_blink = [sub_blink1; sub_blink2];
                    blink_matrix = [blink_matrix; trial_sub_blink]; 
                    trial_sub_opacity = [sub_opacity1; sub_opacity2];
                    current_subject_opacity = [current_subject_opacity; trial_sub_opacity];
                                      
                else

                    sub_blink1 = cell2mat(eyelinkTable.BlinkStublinkQuadrant(row));
                    sub_blink1(1,:) = []; %Remove left eye

                    sub_blink2 = cell2mat(eyelinkTable.BlinkStublinkCenter(row));
                    sub_blink2(1,:) = []; %Remove left eye 
                    
                    sub_opacity1 = eyelinkTable.FaceQuadrantOpacity(row);
                    sub_opacity2 = eyelinkTable.FaceCenterOpacity(row);

                    trial_sub_blink = [sub_blink1; sub_blink2];
                    blink_matrix = [blink_matrix; trial_sub_blink];
                    trial_sub_opacity = [sub_opacity1; sub_opacity2];
                    current_subject_opacity = [current_subject_opacity; trial_sub_opacity];
                    
                end

            end
            
        %Quadrant relevant or irrelevant
        elseif isequal(current_location, 'Quadrant Relevant') 

            %Loop over eyelink rows
            for row = 1:size(eyelinkTable,1)

                % Relevant first
                if eyelinkTable.RelevantFirst(row) == 1

                    sub_blink1 = cell2mat(eyelinkTable.BlinkStublinkQuadrant(row));
                    sub_blink1(1,:) = []; %Remove left eye

                    sub_blink2 = cell2mat(eyelinkTable.BlinkStublinkCenter(row));
                    sub_blink2(1,:) = []; %Remove left eye 
                    
                    sub_opacity1 = eyelinkTable.FaceQuadrantOpacity(row);
                    sub_opacity2 = eyelinkTable.FaceCenterOpacity(row);
                    
                    trial_sub_blink = [sub_blink1; sub_blink2];
                    blink_matrix = [blink_matrix; trial_sub_blink];
                    trial_sub_opacity = [sub_opacity1; sub_opacity2];
                    current_subject_opacity = [current_subject_opacity; trial_sub_opacity];
                    
                else

                    sub_blink1 = cell2mat(eyelinkTable.BlinkStublinkCenter(row));
                    sub_blink1(1,:) = []; %Remove left eye

                    sub_blink2 = cell2mat(eyelinkTable.BlinkStublinkQuadrant(row));
                    sub_blink2(1,:) = []; %Remove left eye 
                    
                    sub_opacity1 = eyelinkTable.FaceCenterOpacity(row);
                    sub_opacity2 = eyelinkTable.FaceQuadrantOpacity(row);
                    
                    trial_sub_blink = [sub_blink1; sub_blink2];
                    blink_matrix = [blink_matrix; trial_sub_blink];
                    trial_sub_opacity = [sub_opacity1; sub_opacity2];
                    current_subject_opacity = [current_subject_opacity; trial_sub_opacity];
                    
                end
                
            end

        end 
                   
        %Extract Run Number
        trial_run_numbers = eyelinkTable.RunNumber;
        
        %% Load Irrelevant Face Event Times
    
        %Define event directory
        if isequal(run_location,'l')
            
        elseif isequal(run_location,'s')
        
            %MRI event directory
            event_dir = ['/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI/',current_ID,'/Perception Task/',current_location,'/MRI Session/MRI Analysis/Event Times'];       

            %Load event times
            cd(event_dir)
            load MRI_jitter_event_times.mat
            
        end

        %% Find Subjects Total PP and PnP Contributions
            
        %Find trials with blinks during blink_rejection_window for current run number
        all_blink_at_stim_trials = find(sum(blink_matrix(:,blink_reject_window),2) > 0);                       

        %Remove trials with blinks during rejection window from BOLD, times, and scores matrices
        subject_eyelink_idx_copy = subject_eyelink_idx;
        subject_eyelink_idx_copy(all_blink_at_stim_trials) = [];
        current_subject_opacity_copy = current_subject_opacity;
        current_subject_opacity_copy(all_blink_at_stim_trials) = [];        
        
        %Find PP and PnP Trials
        PP_idx = subject_eyelink_idx_copy == 1 & current_subject_opacity_copy > 0 & current_subject_opacity_copy < 1;
        PnP_idx = subject_eyelink_idx_copy == 0 & current_subject_opacity_copy > 0 & current_subject_opacity_copy < 1;
        
        all_PP_num = sum(PP_idx);
        all_PnP_num = sum(PnP_idx);       
        
        %Keep PP/PnP number
        keep_PP_num = round(PP_ratio*all_PP_num);
        keep_PnP_num = round(PnP_ratio*all_PnP_num);
        
        %All trials sum
        all_keep_PP_num = all_keep_PP_num + keep_PP_num;
        all_keep_PnP_num = all_keep_PnP_num + keep_PnP_num;   
        
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
                    
                %Run directory
                BOLD_PC_run_data = fullfile('/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI',current_ID,'Perception Task',current_location,...
                    'MRI Session/MRI Analysis/Extracted voxel data',['Run_',num2str(current_run)],'axialSlices/Relevant Stimuli');

                %Session directory 
                BOLD_PC_session_data = fullfile('/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI',current_ID,'Perception Task',current_location,...
                    'MRI Session/MRI Analysis/Extracted voxel data/Session percent change/Relevant Stimuli');
                                    
            end            
            
            %Load BOLD
            load(fullfile(BOLD_PC_run_data,'extracted_signal_PC_jitter_norm_binned_cuts_filtered.mat'))

            %Cut eyelink index and blink matrix to current run trials
            if current_run == 1
                
                extended_trial_idx = [min(find(trial_run_numbers == current_run)):max(find(trial_run_numbers == current_run))*2]';
               
            elseif current_run == 4 && isequal(current_ID,'620') && isequal(current_location, 'Center Relevant') || ...
                current_run == 4 && isequal(current_ID,'622') && isequal(current_location, 'Quadrant Relevant')
            
                extended_trial_idx = [max(find(trial_run_numbers == 2))*2+1:max(find(trial_run_numbers == current_run))*2]';
           
            elseif current_run == 3 && isequal(current_ID,'621') && isequal(current_location, 'Center Relevant') || ...
                    current_run == 3 && isequal(current_ID,'607') && isequal(current_location, 'Quadrant Relevant')
                
                extended_trial_idx = [max(find(trial_run_numbers == 1))*2+1:max(find(trial_run_numbers == current_run))*2]';
                
            elseif current_run == 5 && isequal(current_ID,'607') && isequal(current_location, 'Quadrant Relevant') || ...
                current_run == 5 && isequal(current_ID,'647') && isequal(current_location, 'Quadrant Relevant') || ...
                current_run == 5 && isequal(current_ID,'650') && isequal(current_location, 'Quadrant Relevant') 
            
                extended_trial_idx = [max(find(trial_run_numbers == 3))*2+1:max(find(trial_run_numbers == current_run))*2]';
                              
            else
               
                extended_trial_idx = [max(find(trial_run_numbers == current_run-1))*2+1:max(find(trial_run_numbers == current_run))*2]';
               
            end
            
            %Cut subject confidence score and blink matrix
            run_eyelink_idx = subject_eyelink_idx(extended_trial_idx);
            run_blink_matrix = blink_matrix(extended_trial_idx,:);
            run_opacity = current_subject_opacity(extended_trial_idx,:);

            %Find trials with blinks during blink_rejection_window for current run number
            blink_at_stim_trials = find(sum(run_blink_matrix(:,blink_reject_window),2) > 0);                       
            
            %Remove trials with blinks during rejection window from BOLD, times, and scores matrices
            run_eyelink_idx(blink_at_stim_trials) = [];
            run_opacity(blink_at_stim_trials) = [];
            eval(['jitter_times_run_',num2str(current_run),'(blink_at_stim_trials) = [];'])
            times_img_3D_PC_cut(:,blink_at_stim_trials) = [];
            img_3D_PC_cut(:,:,blink_at_stim_trials) = [];
            
            %Confirm appropriate index and data trial numbers
            if not(isequal(size(run_eyelink_idx,1), size(run_opacity,1), size(img_3D_PC_cut,3), size(times_img_3D_PC_cut,2),...
                    size(eval(['jitter_times_run_',num2str(current_run)]),1)))

                error([current_ID,' -  Number of trials mismatch with BOLD and index'])
                %warning([current_ID,' -  Number of trials mismatch with BOLD and index'])
                continue

            end
            
            %Find PP and PnP Trials
            PP_idx = run_eyelink_idx == 1 & run_opacity > 0 & run_opacity < 1;
            PnP_idx = run_eyelink_idx == 0 & run_opacity > 0 & run_opacity < 1;
            
            %Find predicted perceived and not perceived trials for
            %specified SVM score threshold
            if [size(PP_BOLD,3) + sum(PP_idx)] + [size(PnP_BOLD,3) + sum(PnP_idx)] > 160%185
                
                disp(['Skipping run because too much data'])
                           
                %Clear BOLD data variable
                clearvars img_3D_PC_cut
                
                continue
             
            end
            
            %Check if too many PP trials
            if [size(PP_BOLD,3) + sum(PP_idx)] > 80%100
                
                disp('PP data too large - not adding data for this run')
                
            else
            
                PP_BOLD = cat(3, PP_BOLD, img_3D_PC_cut(:,:,PP_idx));
                PP_BOLD_time = cat(2, PP_BOLD_time, times_img_3D_PC_cut(:,PP_idx));                              
                eval(['PP_event_times_run_',num2str(current_run),' = jitter_times_run_',num2str(current_run),'(PP_idx);'])

            end
            
            %Check if too many PnP trials
            if [size(PnP_BOLD,3) + sum(PnP_idx)] > 80%100
                
                disp('PnP data too large - not adding data for this run')
                
            else
            
                PnP_BOLD = cat(3, PnP_BOLD, img_3D_PC_cut(:,:,PnP_idx));   
                PnP_BOLD_time = cat(2, PnP_BOLD_time, times_img_3D_PC_cut(:,PnP_idx));
                eval(['PnP_event_times_run_',num2str(current_run),' = jitter_times_run_',num2str(current_run),'(PnP_idx);'])

            end
            
            %Clear BOLD data variable
            clearvars img_3D_PC_cut
            
        end    
        
        %Match sample size
        if match_sample_size == 1
            
            %Reset variables
            PP_unfilled_num = 0;
            PnP_unfilled_num = 0;
            
            %Check if the size of PP_BOLD is greater than keep_PP_num
            if size(PP_BOLD,3) > keep_PP_num
                
                %Randomly select subset of trials
                rand_PP = randsample(size(PP_BOLD,3),size(PP_BOLD,3)-keep_PP_num);

                %Cut PP_BOLD to only those randomly keep trials
                PP_BOLD(:,:,rand_PP) = [];
                PP_BOLD_time(:,rand_PP) = [];
                
                %Check sizes
                if not(isequal(size(PP_BOLD,3),size(PP_BOLD,3),size(PP_BOLD_time,2)))
                   
                    error('Size of PP trials is off after random selection!')
                    
                end
                
            % Check if not enough trials
            elseif size(PP_BOLD,3) < keep_PP_num
                
                PP_unfilled_num = [keep_PP_num - size(PP_BOLD,3)];
                warning('Not enough PP trials to keep')
                
            end
            
            %Check if the size of PnP_BOLD is greater than keep_PnP_num
            if size(PnP_BOLD,3) > keep_PnP_num
                
                %Randomly select subset of trials
                rand_PnP = randsample(size(PnP_BOLD,3),size(PnP_BOLD,3)-keep_PnP_num);
                
                %Cut PP_BOLD to only those randomly keep trials
                PnP_BOLD(:,:,rand_PnP) = [];
                PnP_BOLD_time(:,rand_PnP) = [];
                
                %Check sizes
                if not(isequal(size(PnP_BOLD,3),size(PnP_BOLD,3),size(PnP_BOLD_time,2)))
                   
                   error('Size of PnP trials is off after random selection!')
                    
                end
                
            % Check if not enough trials
            elseif size(PnP_BOLD,3) < keep_PnP_num
                
               PnP_unfilled_num = [keep_PnP_num - size(PnP_BOLD,3)];
               warning('Not enough PnP trials to keep')
          
            end
                       
        end
                
        toc
        
        %Save the PP PnP event times
        cd(event_dir) 
%        save([eyelink_type,'_jitter_event_times.mat'],'PP_event*','PnP_event*')
   
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
               
               %Enter zeros trials for PP or PnP count
               if strcmp(current_data_type,'PP_BOLD')
                   
                   PP_trial_num = 0;
               
               elseif strcmp(current_data_type,'PnP_BOLD')
                   
                   PnP_trial_num = 0;
                   
               end
                              
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
        
        if isequal(match_sample_size,1)
        
            save(['extracted_signal_PC_',eyelink_type,'_jitter_matched.mat'], ...
               'PP_BOLD_subject_img_3D_PC_cut', 'PnP_BOLD_subject_img_3D_PC_cut', 'PP_trial_num', 'PnP_trial_num',...
               'PP_unfilled_num', 'PnP_unfilled_num','-v7.3')
           
        else
        
            save(['extracted_signal_PC_',eyelink_type,'_jitter.mat'], ...
               'PP_BOLD_subject_img_3D_PC_cut', 'PnP_BOLD_subject_img_3D_PC_cut', 'PP_trial_num', 'PnP_trial_num', '-v7.3')

        end
    
        toc
        
        %Clear variables before continuing
        clearvars PP_BOLD* PnP_BOLD* subject_img* PP_trial_num PnP_trial_num
        %}
    end  
    
end
