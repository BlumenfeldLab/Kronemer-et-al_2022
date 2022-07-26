%% Create Group Data of EEG Voltage Data - Report Paradigm CP and CnP Epochs

%This code will gather all the CP and CNP trials across subjects, average
%within subjects, and then create a group matrix of channel x time x
%subjects. The code also rejects subjects by poor behavioral and rejects
%trials with a blink event during the stimulus time. 

%Written by: Sharif I. Kronemer
%Date: 5/19/2021

clear

%% Run Location

%Select run location
prompt_1 = 'Running code local or server [l, s]: ';
run_location = input(prompt_1,'s');

%Location relevant
prompt_3 = 'Report post-trial times [1s, 15s, combine]: ';
post_trial_times = input(prompt_3, 's');

%Trials post-stim times
if isequal(post_trial_times,'1s')
    
    poststim_delay = {'1s'};
    delay_save_name = '1s';

elseif isequal(post_trial_times,'15s')
    
    poststim_delay = {'15s'};
    delay_save_name = '15s';
    
elseif isequal(post_trial_times,'combine')
    
    poststim_delay = {'1s','15s'};
    delay_save_name = '1_15s';
    
end

%Blink rejection window (-200 +500)
blink_reject_window = [5801:6500]; %6001:6050

%% Local directories
if isequal(run_location, 's')
    
    %Add group directory to path
    addpath('/mnt/Data8/HNCT sEEG Study/sEEG Analysis/Group Data/Voltage')
    addpath('/mnt/Data6/HNCT sEEG Study/sEEG Analysis/Analysis Code/Behavioral Analysis')
    
    %Subject folder directory
    Report_subject_folder = '/mnt/Data6/HNCT sEEG Study/sEEG Analysis/Subject Analysis';
    
    %Save directory
    save_dir = '/mnt/Data6/HNCT sEEG Study/sEEG Analysis/Group Data/Voltage';
        
elseif isequal(run_location, 'l')
       
end

%Make save directory
mkdir(save_dir)
    
%% Group Variables and Subject Lists

%Find all subject folders
Report_subject_list = dir(Report_subject_folder);
Report_subject_list = {Report_subject_list.name}';

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
for sub = 1:length(Report_subject_list) 
    
    %Select ID
    ID = Report_subject_list{sub};
    
    %Initalize subject variables
    subject_CP_data = [];
    subject_CnP_data = [];
    
    %Skip rejected subjects (Note: 238 - sleeping during study; 367 - EEG
    %data deleted)
    if isequal(ID,'238RC') || isequal(ID,'367NG')
       
       disp(['Skipping ',ID])
       continue
        
    end
    
    %Loop over post stim delay
    for loc = 1:length(poststim_delay)
        
        %Select relevant location
        delay = poststim_delay{loc};
        
        %Rejecting subjects by behavioral performance
        bad_subject_idx = report_subject_rejection_by_behavior_EEG(ID);
        
        %Check if bad subject by behavior
        if isequal(bad_subject_idx, 1)
            
            disp(['Bad behavior - Skipping ', num2str(ID)])
            continue
            
        end

        %Subject percent change data path 
        if isequal(run_location, 'l')
        
        elseif isequal(run_location, 's')
            
            %EEG epochs path
            rootpath = fullfile('/mnt/Data6/HNCT sEEG Study/sEEG Analysis/Subject Analysis/',ID,'/EEG Analysis/Perception Task/Preprocessed Data/Voltage Epochs');
            
            %Bad EEG epochs path
            bad_trial_path = fullfile('/mnt/Data6/HNCT sEEG Study/sEEG Analysis/Subject Analysis/',ID,'/EEG Analysis/Perception Task/Event Times and Identifiers');
            
            %Eyelink path
            eyelink_path = fullfile('/mnt/Data6/HNCT sEEG Study/sEEG Analysis/Subject Analysis/',ID,'/Pupillometry Analysis/Perception Task/Eyelink Table');
            
        end

        %Check if the data mat file exists
        if ~exist(rootpath)
            
            disp(['Cannot find data/directory - Skipping ', num2str(ID),' ',delay])
            continue
            
        else
            
            disp(['**Adding No Report Subject ',num2str(ID),' Post-Stim ', delay,'**'])

        end
        
        %% Load Subject EEG Data - CP and CNP
        
        cd(rootpath) 
        
        %Load CP EEG data
        load(['CP_',delay,'_4s_epochs.mat'])
        CP_epochs = eval(['CP_',delay,'_epochs']);

        %Load CnP EEG data
        load(['CnP_',delay,'_4s_epochs.mat'])
        CnP_epochs = eval(['CnP_',delay,'_epochs']);
            
        %% Exclude trials with blinks at time of stimulus onset
        
        %Load eyelink table
        cd(eyelink_path)
        load('eyelinkTable.mat')

        eyelinkTable(eyelinkTable.FaceQuadrantOpacity == 0,:) = []; %Blank trials
        eyelinkTable(eyelinkTable.FaceQuadrantOpacity == 1,:) = []; %Opaque trials
        eyelinkTable(eyelinkTable.FaceQuadrantOpacity > 1,:) = []; %For 618 & 631 MRI that has a run of opacity greater that 1

       if isequal(delay,'1s')

           eyelinkTable = eyelinkTable(eyelinkTable.PostStimDuration == 1,:);

       elseif isequal(delay,'15s')

           eyelinkTable = eyelinkTable(eyelinkTable.PostStimDuration == 15,:);

       end
                             
       %Subject 366MA is missing Run 1 and 2 EEG files
        if isequal(ID,'366MA')
            
            eyelinkTable(eyelinkTable.RunNumber == 1,:) = []; 
            eyelinkTable(eyelinkTable.RunNumber == 2,:) = []; 
           
        end

        %Create CP and CnP tables
        CP_table = eyelinkTable(eyelinkTable.ConditionCode == 1,:);
        CnP_table = eyelinkTable(eyelinkTable.ConditionCode == 2,:);
        
        %Create blink matrix - Note: all stimuli in this paradigm are
        %quadrant, therefore, it is not necessary to specify center vs
        %quadrant blink as for the no-report paradigm.
        
        %Convert to matrix and remove left eye
        CP_blink_matrix = cell2mat(CP_table.BlinkStublinkQuadrant);
        CP_blink_matrix(1:2:size(CP_table,1)*2,:) = [];    

        CnP_blink_matrix = cell2mat(CnP_table.BlinkStublinkQuadrant);
        CnP_blink_matrix(1:2:size(CnP_table,1)*2,:) = []; 

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
        
        %Remove Run 3 and 4 357HW because the Eyelink edf file for run 3 and 4 is corrupted so there will be a mismatch in EEG
        %and eyelink trial numbers. 
        if isequal(ID,'357HW') && isequal(delay,'1s')
            
            %Remove run 3 and 4
            CP_epochs(:,:,12:26) = [];
            CnP_epochs(:,:,10:20) = [];
            
        elseif isequal(ID,'357HW') && isequal(delay,'15s')  
            
            %Remove run 3 and 4
            CP_epochs(:,:,15:26) = [];
            CnP_epochs(:,:,7:14) = [];
            
        end
                
        %Confirm the number of EEG trials equals Eyelink trials
        if not(isequal(size(CP_epochs,3), size(CP_blink_matrix,1))) 

            error('Number of CP EEG and Eyelink trials mistmatch')
            
        elseif not(isequal(size(CnP_epochs,3), size(CnP_blink_matrix,1))) 

            error('Number of CnP EEG and Eyelink trials mistmatch')
    
        end
        
        %% Trial rejections by EEG and Eyelink parameters
        
        %Load bad trial index
        cd(bad_trial_path)
        load('Bad_trials_index.mat')    
        
        %All trials to reject between EEG and Eyelink rejections
        CP_reject_trials = union(CP_blink_at_stim_trials, eval(['CP_',delay,'_bad_epochs_idx']));
        CnP_reject_trials = union(CnP_blink_at_stim_trials, eval(['CnP_',delay,'_bad_epochs_idx']));
        
        %Remove rejected trials out of range for 357 that had run 3 and 4
        %removed
        if isequal(ID,'357HW') 
            
           %Remove reject trials out of range
           CP_reject_trials((CP_reject_trials > size(CP_epochs,3))) = [];
           CnP_reject_trials((CnP_reject_trials > size(CnP_epochs,3))) = [];
            
        end
                
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
%         CP_epochs(:,:,eval(['CP_',delay,'_bad_epochs_idx'])) = [];
%         CnP_epochs(:,:,eval(['CnP_',delay,'_bad_epochs_idx'])) = [];
        
        %Count the number of bad EEG trials
        CP_reject_trial_num = CP_reject_trial_num + length(CP_reject_trials); %length(eval(['CP_',delay,'_bad_epochs_idx']));
        CnP_reject_trial_num = CnP_reject_trial_num + length(CnP_reject_trials); %length(eval(['CnP_',delay,'_bad_epochs_idx']));
       
        %% Add subject epochs data

        % Create subject matrix of channels x time x trials
        subject_CP_data = cat(3, subject_CP_data, CP_epochs); 
        subject_CnP_data = cat(3, subject_CnP_data, CnP_epochs); 
        
    end
    
    %% Add subject data to group data matrix

    %Add subject CnP data to group 
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

save(['Group_EEG_voltage_CP_CnP_',delay_save_name,'_data.mat'],'CP_epochs_subjects_list','CnP_epochs_subjects_list', ...
    'group_CP_EEG_voltage_data','group_CnP_EEG_voltage_data','CP_reject_trial_num','CnP_reject_trial_num',...
    'CP_total_trial_num','CnP_total_trial_num','-v7.3')
