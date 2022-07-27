%% Test if the irrelevant trial counts are accurate

%The purpose of this code is to confirm that the event times and number of
%fMRI epochs are matching after correcting for rejections by EyeLink,
%behavior, and movement. 

%Written by: Sharif I. Kronemer
%Date: 9/24/2021

%Define excess motion run matrix
motion_dir = 'S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Group Analysis\MRI Analysis\Movement Artifact';

%Define the irrelevant location 
location_list = {'Quadrant','Center'};

%Find all the ID names
ID_list = dir('Y:\HNCT No Report Paradigm\Subject Analysis MRI');
ID_list = {ID_list.name}';
ID_list(1:2) = [];

%Loop over the ID lists
for sub = 1:length(ID_list)
    
    %Select the current ID
    ID = ID_list{sub};
    
    %Loop over the irrelevant locations
    for loc = 1:2
        
        %Initialize matrix
        PP_trial_num = [];
        PnP_trial_num = [];
        PP_combine = [];
        PnP_combine = [];
        
        %Select current location
        location = location_list{loc};
        
        %Load motion artifact vector
        if isequal(loc, 1)
            
            %Center relevant/quadrant irrelevant
            load(fullfile(motion_dir,'MRI_move_artifact_center_rel_2mm_1deg_thresholds.mat'));
            
        else
            
            %Quadrant relevant/center irrelevant
            load(fullfile(motion_dir,'MRI_move_artifact_quadrant_rel_2mm_1deg_thresholds.mat'));
            
        end
        
        %Define the subject-level directory
        sub_dir = ['Y:\HNCT No Report Paradigm\Subject Analysis MRI\',ID,'\Perception Task\', location,' Irrelevant\MRI Analysis'];
        event_dir = fullfile(sub_dir,'Event Times');
        fMRI_dir = fullfile(sub_dir,'Extracted voxel data\Session percent change\Irrelevant Stimuli');
        
        %Load event data        
        try
           
            cd(event_dir)
            load('PP_PnP_irrelevant_face_event_times_confidence_score_0.75.mat')
            
        catch
            
            disp([ID, ' cannot open event variable'])
            continue
            
        end
        
        %Load fMRI data       
        try
            
            cd(fMRI_dir)
            load('extracted_signal_PC_predicted_threshold_opacity_score_0.75.mat','PP_trial_num','PnP_trial_num')
        
        catch
            
            disp([ID, ' cannot open fMRI variable'])
            continue
            
        end
        
        %% Load Subject Eyelink Data
        
        %Quadrant Relevant and Center Irrelevant    
        if isequal(location, 'Center')

            %Eyelink path
            eyelink_path = fullfile('Y:\HNCT No Report Paradigm\Subject Analysis MRI',ID,'Perception Task\Quadrant Relevant\MRI Session\Pupillometry Analysis');

            %Center Relevant and Quadrant Irrelevant
        elseif isequal(location, 'Quadrant')

            %Eyelink path
            eyelink_path = fullfile('Y:\HNCT No Report Paradigm\Subject Analysis MRI',ID,'Perception Task\Center Relevant\MRI Session\Pupillometry Analysis');

        end

        %Load eyelink table
        cd(eyelink_path)
        load('eyelinkTable.mat')

        %Center relevant or irrelevant
        if isequal(location, 'Center')

            %Keep only threshold trials
            eyelinkTable(eyelinkTable.FaceCenterOpacity == 0,:) = []; %Blank trials
            eyelinkTable(eyelinkTable.FaceCenterOpacity == 1,:) = []; %Opaque trials
            eyelinkTable(eyelinkTable.FaceCenterOpacity > 1,:) = []; %For 618 & 631 MRI that has a run of opacity greater that 1

            %Create blink matrix
            
            %Convert to matrix and remove left eye
            blink_matrix = cell2mat(eyelinkTable.BlinkStublinkCenter);
            blink_matrix(1:2:size(eyelinkTable,1)*2,:) = [];   
            
        %Quadrant relevant or irrelevant
        elseif isequal(location, 'Quadrant') 
            
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
        sub_runs = unique(trial_run_numbers); 
                
        %Combine run events
        %Loop over runs
        for run = 1:length(sub_runs)

            current_run = sub_runs(run);

            %Skip bad sessions with excessive movement 
             if  ismember([ID,' Run_',num2str(current_run)],num_bad_runs_by_mm_deg)

                 continue

             end

            %Combine the event times across runs
            PP_combine = [PP_combine; eval(['PP_event_times_run_',num2str(current_run)])]; 
            PnP_combine = [PnP_combine; eval(['PnP_event_times_run_',num2str(current_run)])];

        end
                 
        %Check if trial count matches
        if not(isequal(length(PP_combine),PP_trial_num)) || not(isequal(length(PnP_combine),PnP_trial_num))
            
            disp([ID, ' is wrong - please check!'])
            
        else
            
            disp([ID,' is confirmed!'])
            
        end
        
        %Clear current subject/location variables
        clearvars PP* PnP*
        
    end

end
