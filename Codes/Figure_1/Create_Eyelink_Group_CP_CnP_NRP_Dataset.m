%% Create Group Data of Eyelink Data - No Report CP and CnP Epochs

%This code will gather all the CP and CNP trials across subjects, average
%within subjects, and then create a group matrix of data type x time x
%subjects. The code also rejects subjects by poor behavioral and rejects
%trials with a blink event during the stimulus time. 

%Written by: Sharif I. Kronemer
%Date: 5/20/2021

clear

%% Run Location

%Select run location
prompt_1 = 'Running code local or server [l, s]: ';
run_location = input(prompt_1,'s');

%Location relevant
prompt_2 = 'No-report relevant location [c, q, combine]: ';
relevant_location = input(prompt_2, 's');

%Opacity
prompt_3 = 'Opacity [threshold, opaque]: ';
opacity = input(prompt_3,'s');

%Center Relevant
if isequal(relevant_location, 'c')
    
    relevant_location = 'Center Relevant';
    relevance_cell = {'Center Relevant'};
    rel_save_name = 'cent';

%Quadrant Relevant
elseif isequal(relevant_location, 'q')
    
    relevant_location = 'Quadrant Relevant';
    relevance_cell = {'Quadrant Relevant'};
    rel_save_name = 'quad';
    
%Combine Center and Quadrant Relevant   
elseif isequal(relevant_location, 'combine')
    
    relevant_location = 'Center and Quadrant Relevant';
	relevance_cell = {'Center Relevant','Quadrant Relevant'};
    rel_save_name = 'cent_quad';
    
end

%Imaging modality
prompt_4 = 'Modality type [EEG, MRI, both]: ';
modality = input(prompt_4,'s');

%EEG
if isequal(modality,'EEG')
   
    modality_name = 'EEG';
    modality_cell = {'EEG'};
    
%MRI
elseif isequal(modality,'MRI')
    
    modality_name = 'MRI';
    modality_cell = {'MRI'};
    
%EEG and MRI
elseif isequal(modality,'both')
    
    modality_name = 'EEG_MRI';
    modality_cell = {'EEG','MRI'};
    
end

%Remove trials with blink at stimulus onset (50ms time window)
prompt_5 = 'Remove trials with blink at stim [y,n]: ';
remove_blink_trials = input(prompt_5,'s');

%% Local directories
if isequal(run_location, 's')
    
    %Add group directory to path
    addpath('/mnt/Data16/HNCT_AuditoryNRP_Study/HNCT No Report Paradigm/Analysis//Subject Data/Eyelink mat files')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/Behavioral Analysis')
    
    %Subject folder directory
    noreport_beh_subject_folder = '/mnt/Data16/HNCT No Report Paradigm/Subject Analysis EEG';  
    noreport_MRI_subject_folder = '/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI';  
    
    %Save directory
    save_dir = '/mnt/Data16/HNCT_AuditoryNRP_Study/HNCT No Report Paradigm/Analysis/Group Data/EyeLink';
        
elseif isequal(run_location, 'l')
       
end

%Make save directory
mkdir(save_dir)
    
%% Group Variables and Subject Lists

%Find all subject folders
if isequal(modality, 'EEG')
    
    noreport_subject_list = dir(noreport_beh_subject_folder);
    noreport_subject_list = {noreport_subject_list.name}';
    
elseif isequal(modality, 'MRI')
    
    noreport_subject_list = dir(noreport_MRI_subject_folder);
    noreport_subject_list = {noreport_subject_list.name}';
    
elseif isequal(modality, 'both')
    
    %MRI list
    noreport_MRI_subject_list = dir(noreport_MRI_subject_folder);
    noreport_MRI_subject_list = {noreport_MRI_subject_list.name}';
    
    %EEG list
    noreport_EEG_subject_list = dir(noreport_beh_subject_folder);
    noreport_EEG_subject_list = {noreport_EEG_subject_list.name}';
    
    %Take the union of MRI and EEG subject lists
    noreport_subject_list = union(noreport_MRI_subject_list,noreport_EEG_subject_list);
    
end

%Setup subject count
subject_count = 0;

%Initialized variables
group_CP_pupil_data = [];
group_CnP_pupil_data = [];

group_CP_blink_data = [];
group_CnP_blink_data = [];

group_CP_microsac_data = [];
group_CnP_microsac_data = [];

CP_epochs_subjects_list = {};
CnP_epochs_subjects_list = {};

%Reject trial count
CP_reject_trial_num = 0;
CnP_reject_trial_num = 0;

%Total trial count
CP_total_trial_num = 0;
CnP_total_trial_num = 0;

%% NO REPORT PARADIGM SUBJECTS - Aggregate Eyelink over subjects

disp('Aggregate data across no-report subjects')

%Loop over subject
for sub = 1:length(noreport_subject_list) 
    
    %Select ID
    ID = noreport_subject_list{sub};
    
    %Initalize subject variables
    subject_CP_pupil_data = [];
    subject_CnP_pupil_data = [];
    
    subject_CP_blink_data = [];
    subject_CnP_blink_data = [];
    
    subject_CP_microsac_data = [];
    subject_CnP_microsac_data = [];
    
    %Loop over relevant locations
    for loc = 1:length(relevance_cell)
        
        %Select relevant location
        location = relevance_cell{loc};
        
        %Loop over modality type
        for mod = 1:length(modality_cell)
        
            %Select current modality
            current_modality = modality_cell{mod};
            
            %EEG behavioral exclusions 
            if isequal(current_modality,'EEG')
                
                %Rejecting subjects by behavioral performance
                bad_subject_idx = noreport_subject_rejection_by_behavior_EEG(ID, location);

            %MRI behavior exclusions
            elseif isequal(current_modality,'MRI')
                
                %Rejecting subjects by behavioral performance
                bad_subject_idx = noreport_subject_rejection_by_behavior_MRI(ID, location);
                
            end
            
            %Check if bad subject by behavior
            if isequal(bad_subject_idx, 1)

                disp(['Bad behavior - Skipping ', num2str(ID), ' ', current_modality])
                continue

            end

            %Subject percent change data path 
            if isequal(run_location, 'l')

            elseif isequal(run_location, 's')

                %EEG subject directories
                if isequal(current_modality, 'EEG')

                    %Eyelink path
                    eyelink_path = fullfile('/mnt/Data28/HNCT No Report Paradigm/Subject Analysis EEG',ID,'Perception Task',location,'EEG Session/Pupillometry Analysis');
                
                %MRI subject directories
                elseif isequal(current_modality, 'MRI')

                    %Eyelink path
                    eyelink_path = fullfile('/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI',ID,'Perception Task',location,'MRI Session/Pupillometry Analysis');
                
                end
                
            end
        
            %Load eyelink table
            try cd(eyelink_path)

                load('eyelinkTable.mat')   
                disp(['Adding subject ',ID,' ',location,' ',current_modality])

            catch 

                disp(['Skipping subject',ID,' ',location,' ',current_modality])
                continue

            end

            %% Exclude trials with blinks at time of stimulus onset

            %Center relevant
            if isequal(location, 'Center Relevant')

                %Threshold opacity
                if isequal(opacity,'threshold') 

                    eyelinkTable(eyelinkTable.FaceCenterOpacity == 0,:) = []; %Blank trials
                    eyelinkTable(eyelinkTable.FaceCenterOpacity == 1,:) = []; %Opaque trials
                    eyelinkTable(eyelinkTable.FaceCenterOpacity > 1,:) = []; %For 618 & 631 MRI that has a run of opacity greater that 1

                %Opaque opacity
                elseif isequal(opacity,'opaque') 

                    eyelinkTable(eyelinkTable.FaceCenterOpacity == 0,:) = []; %Blank trials
                    eyelinkTable(eyelinkTable.FaceCenterOpacity < 1,:) = []; %Threshold trials

                end

            %Quadrant relevant
            elseif isequal(location, 'Quadrant Relevant')

                %Threshold opacity
                if isequal(opacity,'threshold')

                    eyelinkTable(eyelinkTable.FaceQuadrantOpacity == 0,:) = []; %Blank trials
                    eyelinkTable(eyelinkTable.FaceQuadrantOpacity == 1,:) = []; %Opaque trials
                    eyelinkTable(eyelinkTable.FaceQuadrantOpacity > 1,:) = []; %For 618 & 631 MRI that has a run of opacity greater that 1

                %Opaque opacity
                elseif isequal(opacity,'opaque')

                    eyelinkTable(eyelinkTable.FaceQuadrantOpacity == 0,:) = []; %Blank trials
                    eyelinkTable(eyelinkTable.FaceQuadrantOpacity < 1,:) = []; %Threshold trials

                end

            end 

            %Create CP and CnP tables
            CP_table = eyelinkTable(eyelinkTable.ConditionCode == 1,:);
            CnP_table = eyelinkTable(eyelinkTable.ConditionCode == 2,:);

            %Create blink matrix for relevant condition

            %Center Relevant
            if isequal(location, 'Center Relevant')

                %Interoplated Pupil Data

                %Convert to matrix and remove left eye
                CP_pupil_matrix = cell2mat(CP_table.interpPupilCenter);
                CP_pupil_matrix(1:2:size(CP_table,1)*2,:) = [];    

                CnP_pupil_matrix = cell2mat(CnP_table.interpPupilCenter);
                CnP_pupil_matrix(1:2:size(CnP_table,1)*2,:) = []; 

                %Blink Data

                %Convert to matrix and remove left eye
                CP_blink_matrix = cell2mat(CP_table.BlinkStublinkCenter);
                CP_blink_matrix(1:2:size(CP_table,1)*2,:) = [];    

                CnP_blink_matrix = cell2mat(CnP_table.BlinkStublinkCenter);
                CnP_blink_matrix(1:2:size(CnP_table,1)*2,:) = []; 

                %Microsaccade Data

                %Convert to matrix and remove left eye
                CP_microsac_matrix = cell2mat(CP_table.MicrosaccadeCenter);
                CP_microsac_matrix = CP_microsac_matrix(2:3:size(CP_table,1)*3,:); %Select right eye   

                CnP_microsac_matrix = cell2mat(CnP_table.MicrosaccadeCenter);
                CnP_microsac_matrix = CnP_microsac_matrix(2:3:size(CnP_table,1)*3,:); %Select right eye

            %Quadrant relevant
            elseif isequal(location, 'Quadrant Relevant')

                %Interoplated Pupil Data

                %Convert to matrix and remove left eye
                CP_pupil_matrix = cell2mat(CP_table.interpPupilQuadrant);
                CP_pupil_matrix(1:2:size(CP_table,1)*2,:) = [];    

                CnP_pupil_matrix = cell2mat(CnP_table.interpPupilQuadrant);
                CnP_pupil_matrix(1:2:size(CnP_table,1)*2,:) = []; 

                %Blink Data

                %Convert to matrix and remove left eye
                CP_blink_matrix = cell2mat(CP_table.BlinkStublinkQuadrant);
                CP_blink_matrix(1:2:size(CP_table,1)*2,:) = [];    

                CnP_blink_matrix = cell2mat(CnP_table.BlinkStublinkQuadrant);
                CnP_blink_matrix(1:2:size(CnP_table,1)*2,:) = []; 

                %Microsaccade Data

                %Convert to matrix and remove left eye
                CP_microsac_matrix = cell2mat(CP_table.MicrosaccadeQuadrant);
                CP_microsac_matrix = CP_microsac_matrix(2:3:size(CP_table,1)*3,:); %Select right eye    

                CnP_microsac_matrix = cell2mat(CnP_table.MicrosaccadeQuadrant);
                CnP_microsac_matrix = CnP_microsac_matrix(2:3:size(CnP_table,1)*3,:); %Select right eye

            end

            %Check that pupil, blink, and microsaccade data match
            if not(size(CP_pupil_matrix,1) == size(CP_blink_matrix,1) && size(CP_pupil_matrix,1) == size(CP_microsac_matrix,1))

               error('Number of CP trials different across eyelink data!')

            elseif not(size(CnP_pupil_matrix,1) == size(CnP_blink_matrix,1) && size(CnP_pupil_matrix,1) == size(CnP_microsac_matrix,1))

               error('Number of CnP trials different across eyelink data!')

            end

            %Find blink at stim trials - CP trials
            if isempty(CP_blink_matrix)

                %Set blink trials matrix to empty
                CP_blink_at_stim_trials = [];

            else

                %Find trials with blink at stimulus
                CP_blink_at_stim_trials = find(sum(CP_blink_matrix(:,6001:6050),2) > 0);

            end

            %Find blink at stim trials - CNP trials
            if isempty(CnP_blink_matrix)

                %Set blink trials matrix to empty
                CnP_blink_at_stim_trials = [];

            else

                %Find trials with blink at stimulus
                CnP_blink_at_stim_trials = find(sum(CnP_blink_matrix(:,6001:6050),2) > 0);

            end

            %% Trial rejections by Eyelink parameters

            %Count all subject trials
            CP_total_trial_num = CP_total_trial_num + size(CP_table,1);
            CnP_total_trial_num = CnP_total_trial_num + size(CnP_table,1);

            %Remove trials with blink at stim onset
            if isequal(remove_blink_trials, 'y')
                
                %Remove trials with blink at sitmulus       
                CP_pupil_matrix(CP_blink_at_stim_trials,:) = [];
                CnP_pupil_matrix(CnP_blink_at_stim_trials,:) = [];

                CP_blink_matrix(CP_blink_at_stim_trials,:) = [];
                CnP_blink_matrix(CnP_blink_at_stim_trials,:) = [];

                CP_microsac_matrix(CP_blink_at_stim_trials,:) = [];
                CnP_microsac_matrix(CnP_blink_at_stim_trials,:) = [];

                %Count the number of trials rejected by blink at stimulus time
                CP_reject_trial_num = CP_reject_trial_num + length(CP_blink_at_stim_trials); 
                CnP_reject_trial_num = CnP_reject_trial_num + length(CnP_blink_at_stim_trials); 

            end
                
            %% Add subject epochs data

            %Create subject matrix trials by time
            subject_CP_pupil_data = cat(1, subject_CP_pupil_data, CP_pupil_matrix);
            subject_CnP_pupil_data = cat(1, subject_CnP_pupil_data, CnP_pupil_matrix);

            subject_CP_blink_data = cat(1, subject_CP_blink_data, CP_blink_matrix);
            subject_CnP_blink_data = cat(1, subject_CnP_blink_data, CnP_blink_matrix);

            subject_CP_microsac_data = cat(1, subject_CP_microsac_data, CP_microsac_matrix);
            subject_CnP_microsac_data = cat(1, subject_CnP_microsac_data, CnP_microsac_matrix);

        end

    end
    
    %% Add subject data to group data matrix

    %Add subject CP data to group - Note: using pupil CP as index because
    %all data types will have the same number trial
    if ~isempty(subject_CP_pupil_data)   

        %Add subject ID to subject list
        CP_epochs_subjects_list = [CP_epochs_subjects_list; ID];

        %Average within subject trials if multiple trials per subject
        if size(subject_CP_pupil_data,1) > 1

            subject_CP_pupil_data = nanmean(subject_CP_pupil_data,1);
            subject_CP_blink_data = nanmean(subject_CP_blink_data,1);
            subject_CP_microsac_data = nanmean(subject_CP_microsac_data,1);

        end

        % Create group matrix of [time x subjects]
        if isempty(group_CP_pupil_data) 

            group_CP_pupil_data(:,1) = subject_CP_pupil_data; 
            group_CP_blink_data(:,1) = subject_CP_blink_data;
            group_CP_microsac_data(:,1) = subject_CP_microsac_data;

        else

            group_CP_pupil_data(:,(end+1)) = subject_CP_pupil_data; 
            group_CP_blink_data(:,(end+1)) = subject_CP_blink_data;
            group_CP_microsac_data(:,(end+1)) = subject_CP_microsac_data;

        end 

    end

    %Add subject CnP data to group - Note: using pupil CnP as index because
    %all data types will have the same number trial
    if ~isempty(subject_CnP_pupil_data)   

        %Add subject ID to subject list
        CnP_epochs_subjects_list = [CnP_epochs_subjects_list; ID];

        %Average within subject trials if multiple trials per subject
        if size(subject_CnP_pupil_data,1) > 1

            subject_CnP_pupil_data = nanmean(subject_CnP_pupil_data,1);
            subject_CnP_blink_data = nanmean(subject_CnP_blink_data,1);
            subject_CnP_microsac_data = nanmean(subject_CnP_microsac_data,1);

        end

        % Create group matrix of [time x subjects]
        if isempty(group_CnP_pupil_data) 

            group_CnP_pupil_data(:,1) = subject_CnP_pupil_data; 
            group_CnP_blink_data(:,1) = subject_CnP_blink_data;
            group_CnP_microsac_data(:,1) = subject_CnP_microsac_data;

        else

            group_CnP_pupil_data(:,(end+1)) = subject_CnP_pupil_data; 
            group_CnP_blink_data(:,(end+1)) = subject_CnP_blink_data;
            group_CnP_microsac_data(:,(end+1)) = subject_CnP_microsac_data;

        end 

    end
     
end

%Save group EyeLink data
cd(save_dir)

if isequal(remove_blink_trials,'y')
    
    save(['Group_eyelink_CP_CnP_',rel_save_name,'_',opacity,'_',modality_name,'_data.mat'],'CP_epochs_subjects_list','CnP_epochs_subjects_list', ...
        'group_CP_pupil_data','group_CnP_pupil_data','group_CP_blink_data','group_CnP_blink_data','group_CP_microsac_data','group_CnP_microsac_data',...
        'CP_reject_trial_num','CnP_reject_trial_num','CP_total_trial_num','CnP_total_trial_num','-v7.3')

else
    
    save(['Group_eyelink_CP_CnP_',rel_save_name,'_',opacity,'_',modality_name,'_no_blink_rejections_data.mat'],'CP_epochs_subjects_list','CnP_epochs_subjects_list', ...
        'group_CP_pupil_data','group_CnP_pupil_data','group_CP_blink_data','group_CnP_blink_data','group_CP_microsac_data','group_CnP_microsac_data',...
        'CP_reject_trial_num','CnP_reject_trial_num','CP_total_trial_num','CnP_total_trial_num','-v7.3')  

end
