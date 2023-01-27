%% Cut BOLD Epochs and Bin 

%Goal of this function is to cut, normalize epochs. This code also manages
%the rejection of runs that exceeded motion thresholds and trials when
%there was a blink at the stimulus time. The blink at stimulus trials are
%saved in the Event Times directory. 

%Written by: Sharif Kronemer
%Written/modified: 5/14/2022

function Percent_change_BOLD_epoch_cutting_function(ID, rootpath, irrel_rootpath, run_num, run_list, relevant_dir, motion_dir, savefolder, irrel_savefolder, beh_data_path, run_location)
     
    %EyeLink table directory
    if isequal(run_location, 's')

        eyelink_dir = ['/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI/',ID,'/Perception Task/',relevant_dir,'/MRI Session/Pupillometry Analysis'];

    elseif isequal(run_location, 'l')

        eyelink_dir = ['Y:\HNCT No Report Paradigm\Subject Analysis MRI\',ID,'\Perception Task\',relevant_dir,'\MRI Session\Pupillometry Analysis'];

    end
    
    %Loop over relevance types - relevant and irrelevant
    for rel = 1:2
        
        %Relevant stimuli
        if isequal(rel,1)
            
            disp('Running - Relevant Stimuli')

            %Relevant trial names
            trial_types = {'blank'};%{'CP_thres', 'CP_opaque', 'CnP_thres'};
            
            %Define session directory
            session_dir = fullfile(rootpath,'Session percent change','Relevant Stimuli');
            
            %Blink trials save dir
            blink_trial_save_dir = savefolder;

            %Make session directory if it does not exist
            if ~exist(session_dir)

                mkdir(session_dir)

            end

        %Irrelevant stimuli
        elseif isequal(rel,2)
            
            disp('Running - Irrelevant Stimuli')

            %Irrelevant trial names
            trial_types = {'blank'};%{'threshold', 'opaque', 'blank'};
            
            %Define session directory
            session_dir = fullfile(irrel_rootpath,'Session percent change','Irrelevant Stimuli');

            %Blink trials save dir
            blink_trial_save_dir = irrel_savefolder;
            
            %Make session directory if it does not exist
            if ~exist(session_dir)

                mkdir(session_dir)

            end

        end

        %Loop over trial types
        for type = 1:size(trial_types,2)

            %Select trial type
            Perception_type = char(trial_types(type)); 

            %Display trial type running
            disp(['Running analysis for ',Perception_type, ' trials'])
                
            for run = 1:run_num 

                %disp(['Cutting epochs - Run ', num2str(run)])
                disp(['Cutting epochs - ',run_list(run).name])

                %Subject specific run data folder
                %datafolder = fullfile(rootpath,['Run_', num2str(run)],'axialSlices'); 
                datafolder = fullfile(rootpath,run_list(run).name,'axialSlices'); 

                %Relevant stimuli
                if isequal(rel,1)

                    %Define run save directory
                    %run_savefolder = fullfile(rootpath,['Run_', num2str(run)],'axialSlices','Relevant Stimuli');
                    run_savefolder = fullfile(rootpath,run_list(run).name,'axialSlices','Relevant Stimuli');

                    %Make session directory if it does not exist
                    if ~exist(run_savefolder)

                        mkdir(run_savefolder)

                    end

                %Irrelevant stimuli
                elseif isequal(rel,2)

                    %Define run save directory
                    %run_savefolder = fullfile(rootpath,['Run_', num2str(run)],'axialSlices','Irrelevant Stimuli');
                    %run_savefolder = fullfile(rootpath,run_list(run).name,'axialSlices','Irrelevant Stimuli');
                    run_savefolder = fullfile(irrel_rootpath,run_list(run).name,'axialSlices','Irrelevant Stimuli');

                    %Make session directory if it does not exist
                    if ~exist(run_savefolder)

                        mkdir(run_savefolder)

                    end

                end

                %% Event Times Epoch Parsing
                % Epochs correspond with the event times of interest and selecting X number
                % of volumes before and after that event. A challenge when defining the
                % epoch times is that the task events are occuring in 'real world' time,
                % while the MRI volume acquisition times are separately defined. Time 0 in
                % real world time = volume 1. To correspond the event of interest with a
                % specific volume - that event time was rounded to the nearest interger and
                % then +1 was applied to correspond to a volume. E.g., an event that
                % occured at 10.75 would be applied to volume # 12. 
                
                % Load the percent change file
                cd(datafolder);
                load('extracted_signal_pc_filtered')

                %Load the behavior event times info
                cd(beh_data_path)

                %Relevant stimuli
                if isequal(rel,1)

                    %Relevant event times
                    load('MRI_CP_CnP_event_times.mat')
                    
                    %Blank event times
                    load('MRI_blank_stim_event_times.mat')

                %Irrelevant stimuli
                elseif isequal(rel, 2)

                    %Irrelevant event times
                    load('MRI_irrelevant_face_event_times.mat');

                end

                pre_face_onset = 20; %Seconds pre task event onset
                post_face_onset = 20; %Seconds post task event onset

                %eval([Perception_type{1} '_Run = ' Perception_type{1} '15s_Face_Onset_' Run_num]); 
                %x = eval(['subject_list_and_event_times(subj_num).' Perception_type{1} '1s_Face_Onset_' Run_num]); % x is the face times for the run
                %face_times = eval([Perception_type,'_face_times_run_',num2str(run)]);
                face_times = eval([Perception_type,'_face_times_',lower(run_list(run).name)]);

                epoch_start = round(face_times - pre_face_onset)+1- 4; %Finds the first volume in epoch //ADDED SOME PADDING
                epoch_end = round(face_times + post_face_onset)+1 +2 ; %Finds the last volume in epoch //ADDED SOME PADDING

                tdim = size(percent_change,1);
                xdim = size(percent_change,2);
                ydim = size(percent_change,3);
                zdim = size(percent_change,4);

                percent_change_2d = reshape(percent_change, tdim, xdim*ydim*zdim); %Creates a 2D matrix [Time x # of voxels]
                img_3D_PC_cut = [];
                times_img_3D_PC_cut = [];

                %Loop over faces
                for face_index = 1:length(face_times) 
                    
                    % added: +2/-4 volumes of "padding" in the cuts before olating
                    % normalize everything to 8 seconds. 

                %% Normalization 

                    current_face = face_times(face_index); %face times

                    %Create NaN matrices to be filled later
                    cut_timecourse = nan(length(epoch_start(face_index):epoch_end(face_index)),xdim*ydim*zdim);
                    %cut_timecourse = nan(47,xdim*ydim*zdim); %Specify 47s
                    

                    times_cut_timecourse = nan(length(epoch_start(face_index):epoch_end(face_index)),1);
                    %times_cut_timecourse = nan(47,1); %Specify 47s 

                    if current_face <= 0 %if the face is before the first volume, skip

                        continue

                 %% NOT SURE IF NECESSARY 7/16/2019
                 %{
                    elseif (current_face-B2face) <= 0 %if b2_pre time is less than 0s take only the values after the first mri acquisition
                        cut_tc = percent_change_2d(1:(epoch_end(face)),:);
                        timevector = (1:(epoch_end(face)))-1-x(face);
                        normalized_timevector = time_normalization(timevector, Q1B1_pre, B1Q2_pre, Q2B2_pre, B2face, 1, Q1B1_post, B1Q2_post, Q2B2_post, ...
                            avg_Q1B1_pre, avg_B1Q2_pre, avg_Q2B2_pre, avg_B2face, 1, avg_Q1B1_post, avg_B1Q2_post, avg_Q2B2_post);
                        %% come back here
                        cut_timecourse((end-size(cut_tc,1)+1):end,:) = cut_tc;
                        times_cut_timecourse((end-size(timevector,1)+1):end,1) = normalized_timevector;
                        disp('Check face times and B2 pre time')
                 %}   

                    %Check if the face epoch start time is greater than the first volume
                    elseif epoch_start(face_index) <= 0 % if the prestim time is after 0 but epoch_start is before 0, then normalize the jitter duration

                        cut_tc = percent_change_2d(1:epoch_end(face_index),:);
                        timevector = (1:(epoch_end(face_index)))-1-face_times(face_index); %transform from volume numbers to seconds

                         %normalized_timevector = time_normalization(timevector, Q1B1_pre, B1Q2_pre, Q2B2_pre, B2face, 1, Q1B1_post, B1Q2_post, Q2B2_post, ...
                %             avg_Q1B1_pre, avg_B1Q2_pre, avg_Q2B2_pre, avg_B2face, 1, avg_Q1B1_post, avg_B1Q2_post, avg_Q2B2_post);
                %        times_cut_timecourse((end-length(timevector)+1):end,1) = normalized_timevector;

                        times_cut_timecourse((end-length(timevector)+1):end,1) = timevector;

                        cut_timecourse((end-size(cut_tc,1)+1):end,:) = cut_tc;

                    %Check if the face epoch end time is less than the total number of volumes    
                    elseif epoch_end(face_index) <= size(percent_change_2d,1) %otherwise, cut out all time points. 
                        
                        cut_tc = percent_change_2d(epoch_start(face_index):epoch_end(face_index),:);
                        timevector = (epoch_start(face_index):epoch_end(face_index))-1-face_times(face_index); %transform from volume numbers to seconds

                %         normalized_timevector = time_normalization(timevector, Q1B1_pre, B1Q2_pre, Q2B2_pre, B2face, 1, Q1B1_post, B1Q2_post, Q2B2_post, ...
                %             avg_Q1B1_pre, avg_B1Q2_pre, avg_Q2B2_pre, avg_B2face, 1, avg_Q1B1_post, avg_B1Q2_post, avg_Q2B2_post);
                %         times_cut_timecourse(:,1) = normalized_timevector;

                        cut_timecourse(:,:) = cut_tc;
                        times_cut_timecourse(:,1) = timevector;

                    end
                    
                    %Fill img_3D_PC_cut with epoch info
                    % img_3D_PC_cut dimensions are (voxels * time * number of face stimuli)
                    if isempty(img_3D_PC_cut)
                    
                        img_3D_PC_cut(:,:,1) = cut_timecourse; %voxel data
                        times_img_3D_PC_cut(:,1) = times_cut_timecourse; %time stamps for voxel data
                    
                    else
                        
                        %Add epoch
                        img_3D_PC_cut(:,:,end+1) = cut_timecourse; %voxel data
                        times_img_3D_PC_cut(:,end+1) = times_cut_timecourse; %time stamps for voxel data
                        
                    end

                    clearvars cut_timecourse

                end
                
                %Save Run Data
                cd(run_savefolder)
                save(['extracted_signal_PC_' Perception_type '_norm_binned_cuts_filtered.mat'],'img_3D_PC_cut','times_img_3D_PC_cut', '-v7.3'); 
    
                %Clear unnecessary variables
                %clearvars -except trial_types beh_data_path datafolder rootpath Perception_type rel run_num session_dir run_list relevant_dir ID motion_dir eyelink_dir blink_trial_save_dir blink_trial_save_dir irrel_rootpath

            end
                
            %% Bin normalized timecourses across runs
             
            %NOTE ON REJECTIONS: This portion of the code will reject runs
            %that exceed motion thresholds and trials when a blink occured
            %during the stimulus
           
            %Load eyelink table to find stimulus time blinks
            load(fullfile(eyelink_dir,'eyelinkTable.mat'))
            
            %Setup the table's CP/CnP condition code value (1/2)
            if isequal(Perception_type,'CP_thres') || isequal(Perception_type,'CP_opaque')
               
                condition_num = 1;
                
            elseif isequal(Perception_type,'CnP_thres') || isequal(Perception_type,'CnP_opaque')
                
                condition_num = 2;
                
            end
            
            %Center relevant
            if isequal(relevant_dir, 'Center Relevant')
                
                %Load motion artifact file
                load(fullfile(motion_dir,'MRI_move_artifact_center_rel_2mm_1deg_thresholds.mat'))
                
                %Relevant stimulus
                if isequal(rel,1)
                
                    %Threshold opacity
                    if isequal(Perception_type,'CP_thres') || isequal(Perception_type,'CnP_thres')

                        eyelinkTable(eyelinkTable.FaceCenterOpacity == 0,:) = []; %Blank trials
                        eyelinkTable(eyelinkTable.FaceCenterOpacity == 1,:) = []; %Opaque trials
                        eyelinkTable(eyelinkTable.FaceCenterOpacity > 1,:) = []; %For 618 & 631 MRI that has a run of opacity greater that 1

                    %Opaque opacity
                    elseif isequal(Perception_type,'CP_opaque') || isequal(Perception_type,'CnP_opaque')

                        eyelinkTable(eyelinkTable.FaceCenterOpacity == 0,:) = []; %Blank trials
                        eyelinkTable(eyelinkTable.FaceCenterOpacity < 1,:) = []; %Threshold trials

                    end
                    
                 %Irrelevant stimulus
                elseif isequal(rel, 2)
                    
                    %Threshold opacity
                    if isequal(Perception_type,'threshold')

                        eyelinkTable(eyelinkTable.FaceQuadrantOpacity == 0,:) = []; %Blank trials
                        eyelinkTable(eyelinkTable.FaceQuadrantOpacity == 1,:) = []; %Opaque trials
                        eyelinkTable(eyelinkTable.FaceQuadrantOpacity > 1,:) = []; %For 618 & 631 MRI that has a run of opacity greater that 1

                    %Opaque opacity
                    elseif isequal(Perception_type,'opaque') 

                        eyelinkTable(eyelinkTable.FaceQuadrantOpacity == 0,:) = []; %Blank trials
                        eyelinkTable(eyelinkTable.FaceQuadrantOpacity < 1,:) = []; %Threshold trials

                    end
                    
                end

            %Quadrant relevant
            elseif isequal(relevant_dir, 'Quadrant Relevant')

                %Load motion artifact file
                load(fullfile(motion_dir,'MRI_move_artifact_quadrant_rel_2mm_1deg_thresholds.mat'))
                  
                %Relevant stimulus
                if isequal(rel,1)
                
                    %Threshold opacity
                    if isequal(Perception_type,'CP_thres') || isequal(Perception_type,'CnP_thres')

                        eyelinkTable(eyelinkTable.FaceQuadrantOpacity == 0,:) = []; %Blank trials
                        eyelinkTable(eyelinkTable.FaceQuadrantOpacity == 1,:) = []; %Opaque trials
                        eyelinkTable(eyelinkTable.FaceQuadrantOpacity > 1,:) = []; %For 618 & 631 MRI that has a run of opacity greater that 1

                    %Opaque opacity
                    elseif isequal(Perception_type,'CP_opaque') || isequal(Perception_type,'CnP_opaque')

                        eyelinkTable(eyelinkTable.FaceQuadrantOpacity == 0,:) = []; %Blank trials
                        eyelinkTable(eyelinkTable.FaceQuadrantOpacity < 1,:) = []; %Threshold trials

                    end
                
                %Irrelevant stimulus
                elseif isequal(rel, 2)
                    
                    %Threshold opacity
                    if isequal(Perception_type,'threshold')

                        eyelinkTable(eyelinkTable.FaceCenterOpacity == 0,:) = []; %Blank trials
                        eyelinkTable(eyelinkTable.FaceCenterOpacity == 1,:) = []; %Opaque trials
                        eyelinkTable(eyelinkTable.FaceCenterOpacity > 1,:) = []; %For 618 & 631 MRI that has a run of opacity greater that 1

                    %Opaque opacity
                    elseif isequal(Perception_type,'opaque') 

                        eyelinkTable(eyelinkTable.FaceCenterOpacity == 0,:) = []; %Blank trials
                        eyelinkTable(eyelinkTable.FaceCenterOpacity < 1,:) = []; %Threshold trials

                    end
                    
                end
                
            end
            
            disp('Bin percent change images per run ...')

            %Generate vector for each bin each vector size is nvoxels*1.
            pre_face_onset = 20; %Seconds pre task event onset
            post_face_onset = 20; %Seconds post task event onset

            %Create empty bins to be filled with BOLD
            for bin = (-pre_face_onset):1:post_face_onset
                
                eval(['bin_time_' num2str(bin+pre_face_onset) ' = [];'])
                
            end

            %For each run, populate each bin with the values.
            for loop_run = 1:run_num
                
                %Subject specific run number
                run = str2double(cell2mat(regexp(run_list(loop_run).name,'\d*','Match')));
                              
                %Reset stim blink trials variables
                blink_at_stim_trials = [];
                
                %Find if run has been flagged with motion artifact and
                %other subject specific rejections
                if any(strcmp([ID,' Run_', num2str(run)], num_bad_runs_by_mm_deg))

                    disp(['Skipping - Run ', num2str(run)])
                    continue
                    
%                 elseif isequal(ID,'607') & isequal(run, 2) & isequal(relevant_dir, 'Quadrant Relevant')
%                     
%                     disp(['Skipping - Run 2'])
%                     continue 
                    
                else
                    
                    disp(['Binning - Run ', num2str(run)])
                    
                end
                
                %% Prepare Eyelink data - Cut out the blink epochs
                
                %Relevant stimulus
                if isequal(rel,1)
   
                    % Blank Trials
                    if isequal(Perception_type,'blank')
                        
                        % Quadrant relevant
                        if isequal(relevant_dir, 'Quadrant Relevant')
                                
                            %Find Blank trials for current run
                            perception_type_table = eyelinkTable(eyelinkTable.RunNumber == run & eyelinkTable.FaceQuadrantOpacity == 0,:);
                        
                        % Center relevant
                        else
                            
                            %Find Blank trials for current run
                            perception_type_table = eyelinkTable(eyelinkTable.RunNumber == run & eyelinkTable.FaceCenterOpacity == 0,:);
                                    
                        end
                    
                    % CP and CnP Trials
                    else
                        
                        %Find CP and CnP trials for current run
                        perception_type_table = eyelinkTable(eyelinkTable.RunNumber == run & eyelinkTable.ConditionCode == condition_num,:);
                    
                    end
                    
                %Irrelevant stimulus
                elseif isequal(rel,2)
                    
                    %Select current run irrelevant trials
                    perception_type_table = eyelinkTable(eyelinkTable.RunNumber == run,:);
                   
                end
                
                %Create blink matrix - Note: need to be careful to select
                %the correct blink data corresponding the relevance and
                %location
                
                %Center relevant
                if isequal(relevant_dir, 'Center Relevant')
                    
                    %Relevant stimulus
                    if isequal(rel,1)

                        %Convert to matrix and remove left eye
                        blink_matrix = cell2mat(perception_type_table.BlinkStublinkCenter);
                        blink_matrix(1:2:size(perception_type_table,1)*2,:) = [];  

                    %Irrelevant stimulus
                    elseif isequal(rel,2)

                        %Convert to matrix and remove left eye
                        blink_matrix = cell2mat(perception_type_table.BlinkStublinkQuadrant);
                        blink_matrix(1:2:size(perception_type_table,1)*2,:) = [];  

                    end
                
                %Quadrant relevant
                elseif isequal(relevant_dir, 'Quadrant Relevant')

                    %Relevant stimulus
                    if isequal(rel,1)

                        %Convert to matrix and remove left eye
                        blink_matrix = cell2mat(perception_type_table.BlinkStublinkQuadrant);
                        blink_matrix(1:2:size(perception_type_table,1)*2,:) = [];  

                    %Irrelevant stimulus
                    elseif isequal(rel,2)

                        %Convert to matrix and remove left eye
                        blink_matrix = cell2mat(perception_type_table.BlinkStublinkCenter);
                        blink_matrix(1:2:size(perception_type_table,1)*2,:) = [];  

                    end
                
                end

                %% Load the PC run trial data
                
                %Relevant stimuli
                if isequal(rel,1)

                    %Define run save directory
                    run_savefolder = fullfile(rootpath,run_list(loop_run).name,'axialSlices','Relevant Stimuli');

                %Irrelevant stimuli
                elseif isequal(rel,2)

                    %Define run save directory
                    run_savefolder = fullfile(irrel_rootpath,run_list(loop_run).name,'axialSlices','Irrelevant Stimuli');

                end
                
                %Load PC cut data (dimensions = time * voxels (XYZ) * faces; time vector dimensions: time * faces)
                cd(run_savefolder)
                load(['extracted_signal_PC_' Perception_type '_norm_binned_cuts_filtered.mat'])

                %Confirm the number of BOLD trials equals Eyelink trials
                if not(isequal(size(times_img_3D_PC_cut,2), size(blink_matrix,1)))
                   
                    error('Number of BOLD and Eyelink trials mistmatch')
                    
                end
                
                %% Load over faces index
                                
                %Find the number of faces and loop over faces
                for face_index = 1:size(times_img_3D_PC_cut,2)
                    
                    %Find if blink event during stimulus (face time =
                    %6001:6050) - Skip to next face if blink at stim time
                    if sum(blink_matrix(face_index,6001:6050)) > 0
                        
                        %Store trials with blink at stim time
                        blink_at_stim_trials = [blink_at_stim_trials, face_index];                       
                        continue

                    end

                    %Find the number to timepoints/volumes and loop over time
                    for vol_time_index = 1:size(times_img_3D_PC_cut,1) 

                        %Find the specific volume time for a particular trial
                        vol_time = times_img_3D_PC_cut(vol_time_index,face_index);

                        %Loop over bins/time
                        for bin = (-pre_face_onset):1:post_face_onset

                            %Find the volume times nearest bin
                            if (vol_time>=bin) && (vol_time<(bin+1)) %new bining method 
                               %(vol_time>(bin-0.5)) && (vol_time<=(bin+0.5))%old bining method 

                                %Store XYZ MR data into specific bin matrix
                                eval(['bin_time_' num2str(bin+pre_face_onset) '(:,end+1) = img_3D_PC_cut(vol_time_index,:,face_index);'])

                            end

                        end

                    end
                    
                end
                
                %Rename stim blink trials
                eval(['stim_blink_trials_',Perception_type,'_run_',num2str(run),' = blink_at_stim_trials;']);
                
            end

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
                        
            %Save PC data
            cd(session_dir)
            save(['extracted_signal_PC_' Perception_type '_norm_binned_cuts_filtered_nostimblinks.mat'],'subject_img_3D_PC_cut','-v7.3')

        end
                    
        %Save index of trials with blink at stim time
        cd(blink_trial_save_dir)
        
        %If Blank Trials
        if isequal(Perception_type,'blank')
            
            save(['Blink_at_blank_stim_trials_index.mat'],'stim_blink_trials*')

        %If Non-Blank Trials
        else
            
            save(['Blink_at_stim_trials_index.mat'],'stim_blink_trials*')
        
        end
    
        clearvars stim_blink_trials*
    
    end
    
end