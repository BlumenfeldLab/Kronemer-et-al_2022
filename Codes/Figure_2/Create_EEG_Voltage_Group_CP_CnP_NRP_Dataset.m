%% Create Group Data of EEG Voltage Data - No Report CP and CnP Epochs

%This code will gather all the CP and CNP trials across subjects, average
%within subjects, and then create a group matrix of channel x time x
%subjects. The code also rejects subjects by poor behavioral and rejects
%trials with a blink event during the stimulus time. 

%Written by: Sharif I. Kronemer
%Date: 5/20/2021
%Modified: 6/8/2021

clear

%% Run Location

%Select run location
% prompt_1 = 'Running code local or server [l, s]: ';
% run_location = input(prompt_1,'s');
run_location = 's';

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

%Blink rejection window in Eyelink timecourse units -6000 to 6000
blink_reject_window = [5801:6500]; %Stimulus time: 6001:6050; Brief window: 5801:6500 (used for main figure); Long window: 5801:8000

%% Local directories
if isequal(run_location, 's')
    
    %Add group directory to path
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Data/EEG')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/Behavioral Analysis')
    
    %Subject folder directory
    Noreport_subject_folder = '/mnt/Data28/HNCT No Report Paradigm/Subject Analysis EEG';  
    
    %Save directory
    save_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Data/EEG/Voltage';
        
elseif isequal(run_location, 'l')
       
end

%Make save directory
mkdir(save_dir)
    
%% Group Variables and Subject Lists

%Find all subject folders
Noreport_subject_list = dir(Noreport_subject_folder);
Noreport_subject_list = {Noreport_subject_list.name}';

%Setup subject count
subject_count = 0;

%Initialized variables
group_CP_EEG_voltage_data = [];
group_CnP_EEG_voltage_data = [];
CP_epochs_subjects_list = {};
CnP_epochs_subjects_list = {};

%Reject trial count
CP_reject_trial_num = 0;
CnP_reject_trial_num = 0;

%Total trial count
CP_total_trial_num = 0;
CnP_total_trial_num = 0;

%% NO REPORT PARADIGM SUBJECTS - Aggregate voltage over subjects

disp('Aggregate data across no-report subjects')

%Loop over subject
for sub = 1:length(Noreport_subject_list) 
    
    %Select ID
    ID = Noreport_subject_list{sub};
    
    %Initalize subject variables
    subject_CP_data = [];
    subject_CnP_data = [];
    
    %Loop over relevant locations
    for loc = 1:length(relevance_cell)
        
        %Select relevant location
        location = relevance_cell{loc};
                
        %Rejecting subjects by behavioral performance
        bad_subject_idx = noreport_subject_rejection_by_behavior_EEG(ID, location);
        
        %Check if bad subject by behavior
        if isequal(bad_subject_idx, 1)
            
            disp(['Bad behavior - Skipping ', num2str(ID)])
            continue
            
        end

        %Subject percent change data path 
        if isequal(run_location, 'l')
        
        elseif isequal(run_location, 's')
            
            %EEG epochs path
            rootpath = fullfile('/mnt/Data28/HNCT No Report Paradigm/Subject Analysis EEG',ID,'Perception Task',location,'EEG Session/EEG Analysis/Preprocessed Data/Voltage Epochs');
            
            %Bad EEG epochs path
            bad_trial_path = fullfile('/mnt/Data28/HNCT No Report Paradigm/Subject Analysis EEG',ID,'Perception Task',location,'EEG Session/EEG Analysis/Event Times and Identifiers');
            
            %Eyelink path
            eyelink_path = fullfile('/mnt/Data28/HNCT No Report Paradigm/Subject Analysis EEG',ID,'Perception Task',location,'EEG Session/Pupillometry Analysis');

        end

        %Check if the data mat file exists
        if ~exist(rootpath)
            
            disp(['Cannot find data/directory - Skipping ', num2str(ID),' ',location])
            continue
            
        else
            
            disp(['**Adding No Report Subject ',num2str(ID),' ', location,'**'])

        end
        
        %% Load Subject EEG Data - CP and CNP
        
        cd(rootpath) 
        
        %Load CP EEG data
        
        %Threshold
        if isequal(opacity,'threshold')
        
            load('CP_threshold_4s_epochs.mat')
            CP_epochs = CP_threshold_epochs;
            
        %Opaque
        elseif isequal(opacity,'opaque')
            
            load('CP_opaque_4s_epochs.mat')
            CP_epochs = CP_opaque_epochs;
            
        end
        
        %Load CnP EEG data
        
        %Threshold
        if isequal(opacity,'threshold')
        
            load('CnP_threshold_4s_epochs.mat')
            CnP_epochs = CnP_threshold_epochs;
            
        %Opaque
        elseif isequal(opacity,'opaque')
            
            load('CnP_opaque_4s_epochs.mat')
            CnP_epochs = CnP_opaque_epochs;
            
        end
        
        %% Exclude trials with blinks at time of stimulus onset
        
        %Load eyelink table
        cd(eyelink_path)
        load('eyelinkTable.mat')

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
            
            %Convert to matrix and remove left eye
            CP_blink_matrix = cell2mat(CP_table.BlinkStublinkCenter);
            CP_blink_matrix(1:2:size(CP_table,1)*2,:) = [];    

            CnP_blink_matrix = cell2mat(CnP_table.BlinkStublinkCenter);
            CnP_blink_matrix(1:2:size(CnP_table,1)*2,:) = []; 
            
        %Quadrant relevant
        elseif isequal(location, 'Quadrant Relevant')
        
            %Convert to matrix and remove left eye
            CP_blink_matrix = cell2mat(CP_table.BlinkStublinkQuadrant);
            CP_blink_matrix(1:2:size(CP_table,1)*2,:) = [];    

            CnP_blink_matrix = cell2mat(CnP_table.BlinkStublinkQuadrant);
            CnP_blink_matrix(1:2:size(CnP_table,1)*2,:) = []; 

        end
        
        %Find blink at stim trials - CP
        if isempty(CP_blink_matrix)
        
            %Set blink trials matrix to empty
            CP_blink_at_stim_trials = [];
            
        else
            
            %Find trials with blink at stimulus
            CP_blink_at_stim_trials = find(sum(CP_blink_matrix(:,blink_reject_window),2) > 0);
            
        end
        
        %Find blink at stim trials - CNP
        if isempty(CnP_blink_matrix)
        
            %Set blink trials matrix to empty
            CnP_blink_at_stim_trials = [];
            
        else
            
            %Find trials with blink at stimulus
            CnP_blink_at_stim_trials = find(sum(CnP_blink_matrix(:,blink_reject_window),2) > 0);
            
        end
        
        %Removed EEG data for Subject 652 Center Relevant Run 3 not available for Eyelink
        if isequal(ID,'652') && isequal(location,'Center Relevant') && isequal(opacity,'threshold')
            
            %Remove run 3 CP and CnP trials
            CP_epochs(:,:,16:21) = [];
            CnP_epochs(:,:,14:23) = [];
            
        elseif isequal(ID,'652') && isequal(location,'Center Relevant') && isequal(opacity,'opaque')
            
            %Remove run 3 CP and CnP trials
            CP_epochs(:,:,4:5) = [];
            CnP_epochs(:,:,2) = [];
            
        end
                   
        %Confirm the number of EEG trials equals Eyelink trials
        if not(isequal(size(CP_epochs,3), size(CP_blink_matrix,1)))

            error('Number of CP EEG and Eyelink trials mistmatch')
            
        elseif not(isequal(size(CnP_epochs,3), size(CnP_blink_matrix,1)))

            error('Number of CnP EEG and Eyelink trials mistmatch')
    
        end
        
        %% Trial rejections by EEG and Eyelink parameters
        
        %Load bad EEG trial index
        cd(bad_trial_path)
        load('Bad_trials_index.mat')
        
        %Subject 652 Center Relevant Run 3 not available for Eyelink -
        %Remove or update the bad trial 
        if isequal(ID,'652') && isequal(location,'Center Relevant') && isequal(opacity,'threshold')
            
            %Replace the bad epochs 35,36 with 25,26 for after run 3 is
            %rejected that exclude 10 CnP trials 
            CnP_threshold_bad_epochs_idx = [2;4;25;26];
            
            %Note: CP bad EEG epochs [11;14] need not be adjusted because they
            %occur before run 3 so that index is maintained for those
            %trials.
            
        end
               
        %All trials to reject between EEG and Eyelink rejections
        CP_reject_trials = union(CP_blink_at_stim_trials, eval(['CP_',opacity,'_bad_epochs_idx']));
        CnP_reject_trials = union(CnP_blink_at_stim_trials, eval(['CnP_',opacity,'_bad_epochs_idx']));
        
        %Count all subject trials
        CP_total_trial_num = CP_total_trial_num + size(CP_epochs,3);
        CnP_total_trial_num = CnP_total_trial_num + size(CnP_epochs,3);
        
        %Remove blink at stim trials
        CP_epochs(:,:,CP_reject_trials) = [];
        CnP_epochs(:,:,CnP_reject_trials) = [];
                
%         %Remove blink at stim trials
%         CP_epochs(:,:,CP_blink_at_stim_trials) = [];
%         CnP_epochs(:,:,CnP_blink_at_stim_trials) = [];
            
%         %Remove bad trials
%         CP_threshold_epochs(:,:,CP_threshold_bad_epochs_idx) = [];
%         CnP_threshold_epochs(:,:,CnP_threshold_bad_epochs_idx) = [];
        
        %Count the number of bad EEG trials
        CP_reject_trial_num = CP_reject_trial_num + length(CP_reject_trials); 
        CnP_reject_trial_num = CnP_reject_trial_num + length(CnP_reject_trials); 
 
        %% Add subject epochs data

        % Create subject trials to subject matrix of channels x time x trials
        subject_CP_data = cat(3, subject_CP_data, CP_epochs); 
        subject_CnP_data = cat(3, subject_CnP_data, CnP_epochs); 
        
    end
    
    %% Add subject data to group data matrix

    %Add subject CP data to group 
    if ~isempty(subject_CP_data)   
        
        %Add subject ID to subject list
        CP_epochs_subjects_list = [CP_epochs_subjects_list; ID];
        
        %Average within subject trials if multiple trials per subject
        if size(subject_CP_data,3) > 1

            subject_CP_data = nanmean(subject_CP_data,3);

        end

        % Create group matrix of [channel x time x subjects]
        if isempty(group_CP_EEG_voltage_data) 

            group_CP_EEG_voltage_data(:,:,1) = subject_CP_data; 

        else

            group_CP_EEG_voltage_data(:,:,(end+1)) = subject_CP_data;

        end 

    end

    %Add subject CnP data to group 
    if ~isempty(subject_CnP_data)
        
        %Add subject ID to subject list
        CnP_epochs_subjects_list = [CnP_epochs_subjects_list; ID];

        %Average within subject if multiple trials per subject        
        if size(subject_CnP_data,3) > 1

            subject_CnP_data = nanmean(subject_CnP_data,3);

        end

        % Create group matrix of [time x voxels x subjects]
        if isempty(group_CnP_EEG_voltage_data) 

            group_CnP_EEG_voltage_data(:,:,1) = subject_CnP_data; 

        else

            group_CnP_EEG_voltage_data(:,:,(end+1)) = subject_CnP_data;

        end

    end
      
end

%Save group EEG data
cd(save_dir)
save(['Group_EEG_voltage_CP_CnP_',rel_save_name,'_',opacity,'_data.mat'],'CP_epochs_subjects_list','CnP_epochs_subjects_list','relevant_location', ...
    'group_CP_EEG_voltage_data','group_CnP_EEG_voltage_data','CP_reject_trial_num','CnP_reject_trial_num',...
    'CP_total_trial_num','CnP_total_trial_num','-v7.3')

%Save name if using the long blink window
% save(['Group_EEG_voltage_CP_CnP_',rel_save_name,'_',opacity,'_data_longblinkwindow.mat'],'CP_epochs_subjects_list','CnP_epochs_subjects_list','relevant_location', ...
%     'group_CP_EEG_voltage_data','group_CnP_EEG_voltage_data','CP_reject_trial_num','CnP_reject_trial_num',...
%     'CP_total_trial_num','CnP_total_trial_num','-v7.3')