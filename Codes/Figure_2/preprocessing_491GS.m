%% This is for getting channel labels after loading from EEGlab
% labels = {};
%
% for i = 1:length(EEG.chanlocs)
%Extract epoch from icEEG data clippings

%Date: 4/16/2019

% DEPENDECIES:
% edfread.m
% find_pulses.m
% read_mixed_csv.m
% verify_psychopy_ttl.m
% face_pulse.mat
% question_pulse.mat
% button_pulse.mat
% Subject's .edf files
% Subject's behavioral data .csv files

clear


prompt_1 = 'Running code local or server [l, s]: ';
run_location = input(prompt_1,'s');

%Add function directory to path
addpath('\\blumenfeld1.med.yale.internal\Data7\RNS Study\Analysis\Analysis Code\Epoch extraction functions')

%% Location of EEG folder/files

eeg_root = '\\blumenfeld1.med.yale.internal\Data7\RNS Study\Subject Raw Data\491GS\Scalp EEG Data\'; 
eeg_folders = {'Calibration and All Runs'};
eeg_files = {'Blumenfeld~ Re_aef83f00-c5e2-43b3-b146-a2d42f93fa7d-archive.lay'};

%% Location of functions
cleanline_dir = '\\blumenfeld1.med.yale.internal\Data7\RNS Study\Analysis\Analysis Code\Preprocessing\eeglab14_0_0b\plugins'; %Cleanline function folder 

%% Parameters
sampling_rate = 256; 
%Voltage threshold (+/- 100 mV)
voltage_threshold = 100; 
%Proportion of high voltage samples
high_volt_prop = 0.2;
prop_bad_channels = 0.25;

%% Location of behavior folder/files

beh_root = {'\\blumenfeld1.med.yale.internal\Data7\RNS Study\Subject Raw Data\491GS\Behavioral Data'};
beh_folders = {'Calibration and All Runs'};
beh_files = {'491GS_Session 1_2019_Jun_13_1334_12.csv'};

%% Define epoch types
epochtypes = {'Face_presentations','Questions','Buttons'};

%Loop over sessions
for session = 1%:length(eeg_folders)
    
    % Location to store epochs
    save_dir = '\\blumenfeld1.med.yale.internal\Data7\RNS Study\Analysis\Subject Analysis\491GS\Scalp EEG Analysis';
    
%     %Create save folder if it does not exist
%     if ~exist(char(fullfile(save_dir, beh_folders(session))))
%         
%         mkdir(char(fullfile(save_dir, beh_folders(session))));
%         
%     end     
%     
%     %Define new save directory
%     save_dir = char(fullfile(save_dir, beh_folders(session)));
    
    %% Load EEG files
    cd(char(fullfile(eeg_root, eeg_folders(session))))
    
    %Utilize the layread function
    disp(['Load EEG data'])
        
    [header record] = layread(char(fullfile(eeg_root, eeg_folders(session), eeg_files(session))));

    %Remove data prior to run phases
    record = record(:,987081:2507176);
    
%{    
    %% Re-reference - Average Reference
    
    %Find mean voltage across EEG channels
    mean_sample = mean(record([1:4,7:22,29:30],:),1);
    
    %Initialize the variable
    avg_record = [];

    %Loop through all the samples in the raw file
    for sample = 1:size(record,2)
        
        %Substract mean value per sample from the value in each channel for that sample
        avg_record(:,sample)= record(:,sample) - mean_sample(:,sample); 

    end
%}
    
    %Select the data from the research channels 
    ttl = record(27,:)-record(26,:); 
    
    %60Hz filter TTL channel
    %Sampling rate
    fs = 256;
    
    %Design filter
    d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',fs);
    
    %Run filter
    ttl = filtfilt(d, ttl); 

        %% Session Level Preprocessing 
    % 0.1Hz filter
    % 60Hz cleanline
    % finds indices for high SD and high freq analyses which are saved in function
    
    EEG = record; 
    ID = '491GS'; 
    cd('\\blumenfeld1.med.yale.internal\data7\RNS Study\Analysis\Analysis Code\Preprocessing\Preprocessing Pipeline Functions')
    [EEG_highpass_cleanline,EEG_blank] = eeg_session_level_processing_RNS(EEG, ID, save_dir, cleanline_dir, run_location, sampling_rate);
        
    %% Run behavioral analysis on the excel data
    
    % Load the relevant file
    disp(['Load behavioral data'])
    [num,text,raw] = xlsread(char(fullfile(beh_root, beh_folders(session), beh_files(session))));
    
    % Remove row that is NaN
    if isnan(cell2mat(raw(2,1)))
       
       raw(2,:) = [];
       
    end
    
    %Find number of rows and columns
    row = length(raw(:,1));
    column = length(raw(1,:));

    %% Add columns to the end of data sheet for calculations
    raw{1,column+1} = 'Perception Accuracy';
    raw{1,column+2} = 'Perception RT';
    raw{1,column+3} = 'Location Accuracy';
    raw{1,column+4} = 'Location RT';

    %% Find the appropriate columns and save to variables
    for i = 1:size(raw,2)
        
        if strcmp('TRIAL TYPE',raw{1,i})
            trialtype = i;
        elseif strcmp('Delay',raw{1,i})
            delay = i;
        elseif strcmp('Perception answer',raw{1,i})
            perceptionanswer = i;
        elseif strcmp('Perception Accuracy',raw{1,i})
            perception = i;
        elseif strcmp('Perception question time',raw{1,i})
            questionp = i;       
        elseif strcmp('Location answer',raw{1,i})
            locationanswer = i;
        elseif strcmp('Location Accuracy',raw{1,i})
            location = i;
        elseif strcmp('Location question time',raw{1,i})
            questionl = i; 
        elseif strcmp('FaceDrawStart',raw{1,i})
            facedraw = i;
        elseif strcmp('BLOCK NUMBER',raw{1,i})
            block = i;
        elseif strcmp('Perception keypress time',raw{1,i})
            buttonp = i;
        elseif strcmp('Location keypress time',raw{1,i})
            buttonl = i;
        elseif strcmp('Perception RT',raw{1,i})
            reactiontimep = i;
        elseif strcmp('Location RT',raw{1,i})
            reactiontimel = i;
        elseif strcmp('Trial start time',raw{1,i})
            trialstart = i;
        elseif strcmp('TRIAL',raw{1,i})
            trial = i;
        elseif strcmp('Face shown',raw{1,i})
            faceshown = i;
        elseif strcmp('Face quadrant',raw{1,i})
            quadrant = i;    
        elseif strcmp('Face opacity',raw{1,i})
            opacity = i;
        end
        
    end

    totalblocks = max(cell2mat(raw(2:end,block)));

    %% Fill in the four columns on perception and location accuracy and reaction times that you created
    
    %Loop over raw variable row
    for j = 2:size(raw,1)
        
        %Use only the rows with actual run trials, omit calibration rows 
        if strcmp(raw{j,trialtype},'NOISE') || strcmp(raw{j,trialtype},'MOVIE') || any(ismember(raw{j,trialtype},[4 2]))

            %Fill in the perception accuracy for that row
            if raw{j,faceshown} == 1 && raw{j,perceptionanswer} == 1
                raw{j,column+1} = 'TP';
            elseif raw{j,faceshown} == 1 && raw{j,perceptionanswer} == 0
                raw{j,column+1} = 'FN';
            elseif raw{j,faceshown} == 0 && raw{j,perceptionanswer} == 1
                raw{j,column+1} = 'FP';
            elseif raw{j,faceshown} == 0 && raw{j,perceptionanswer} == 0
                raw{j,column+1} = 'TN';
            end

            %Fill in perception reaction time for that row
            raw{j,column+2} = raw{j,buttonp} - raw{j,questionp};

            %Fill in location accuracy for that row
            if strcmp(raw{j,perception},'TP') && strcmp(raw{j,quadrant},raw{j,locationanswer}) == 1 
                raw{j,column+3} = 'Confirmed perceived';
            elseif strcmp(raw{j,perception},'FN') && strcmp(raw{j,quadrant},raw{j,locationanswer}) == 0
                raw{j,column+3} = 'Confirmed not perceived';
            elseif strcmp(raw{j,perception},'FN') && strcmp(raw{j,quadrant},raw{j,locationanswer}) == 1
                raw{j,column+3} = 'Correct guess';
            elseif strcmp(raw{j,perception},'TP') && strcmp(raw{j,quadrant},raw{j,locationanswer}) == 0 
                raw{j,column+3} = 'False perception';
            elseif strcmp(raw{j,perception},'TN') 
                raw{j,column+3} = 'True Negative';
            elseif strcmp(raw{j,perception},'FP') 
                raw{j,column+3} = 'False Positive';
            end              

            %Fill in location reaction time for that row
            raw{j,column+4} = raw{j,buttonl} - raw{j,questionl}; 
        end
    end 

    %There is one face presentation per trials
    psychopy_faces = NaN(1,size(raw,1)-1);
    
    %There are 2 questions and 2 button presses per trial
    psychopy_questions = NaN(2,size(raw,1)-1);
    psychopy_buttons = NaN(2,size(raw,1)-1);
    
    for j = 2:size(raw,1)
        
        psychopy_faces(1,j-1) = raw{j,facedraw};
        psychopy_questions(1,j-1) = raw{j,questionp};
        psychopy_questions(2,j-1) = raw{j,questionl};
        psychopy_buttons(1,j-1) = raw{j,buttonp};
        psychopy_buttons(2,j-1) = raw{j,buttonl};
        
    end
        
%% Find the TTL pulses of interest and ID event onset times
    
    %Create time variable 
    time = 0:1/fs:(size(EEG_highpass_cleanline.data, 2) - 1)/fs;
    
    
    %Define width for pulse of interest
    pulse_of_interest = 0.15; % We want 150ms pulses
    
    %Define the range of pulse width 
    min_pulse_length = 0.025;
    max_pulse_length = 0.3;
    
    tol = 0.005;
    
    %Pick a point (voltage value) that cuts through all the pulses
    query_point = 150;
    
    %Run find pulse function
    [onsets, offsets] = find_pulses_CJM(ttl, fs, min_pulse_length, max_pulse_length, query_point, tol);
    
    %Add TTL onset times (specific for 491GS***)
    %onsets  = sort([onsets, 1061949]);
        
    %Add TTL offset times (specific for 491GS***)
    %offsets  = sort([offsets, 1061988]);
    
    %Remove TTL onset times that are not face events (specfic for 491GS***)
    bad_TTL_pulse_idx = find(ismember(onsets,[147534,288634,348766,467075,652105,788282,1202564,1326839,1517930]));
    
    %Remove from onset and offset times
    onsets(bad_TTL_pulse_idx) = [];
    offsets(bad_TTL_pulse_idx) = [];
    
    %Find the duration of the TTL pulse
    durations = time(offsets) - time(onsets);
    pulse_of_interest_onsets = find(abs(durations - pulse_of_interest) <= tol);
        
%     %Real Time - Plot TTL channel and the points marking event onset
%     figure;
%     plot(time, ttl); hold on;
%     %plot(time(onsets), zeros(1, length(onsets)), '*r')
%     
%     plot(time(onsets(pulse_of_interest_onsets)), zeros(1, numel(pulse_of_interest_onsets)),'*r')
%     xlabel('TTL Time')
%     ylabel('TTL signal amplitude (mV)')
%     title('TTL face pulses (in red) during task')
%     
%     %Sample Rate Time - Plot TTL channel and the points marking event onset
%     figure;
%     plot(ttl); hold on;
%     %plot(time(onsets), zeros(1, length(onsets)), '*r')
%     
%     plot(onsets(pulse_of_interest_onsets), zeros(1, numel(pulse_of_interest_onsets)),'*r')
%     xlabel('TTL Time')
%     ylabel('TTL signal amplitude (mV)')
%     title('TTL face pulses (in red) during task')
%     
%     %Sample Rate Time - Plot TTL pulse onset and offset times
%     figure;
%     plot(ttl); hold on;
%     %plot(time(onsets), zeros(1, length(onsets)), '*r')
%     
%     plot(offsets, zeros(1, numel(offsets)),'*r')
%     plot(onsets, zeros(1, numel(onsets)),'*b')
%     xlabel('TTL Time')
%     ylabel('TTL signal amplitude (mV)')
%     title('TTL face pulses (in red) during task')

    %% Extract Epochs
    disp('Cutting epochs')
    
    %Sampling frequency
    ttl_diode_offset = floor(.042 * fs);

    %% Cut short epochs -1 to +1
   
    %Loop over the epoch types
    for type = 1%:length(epochtypes)
        
        event_type = epochtypes{type};
        epochs = [];
        
        %If face epoch type
        if isequal(event_type, 'Face_presentations')
            
           IDX = onsets(pulse_of_interest_onsets);
              
        end    
        
        %Loop over the IDX variable       
        for i = 1:length(IDX)
            
            %Define the start and stop time of the epoch
            attempted_start_idx = IDX(i)+ttl_diode_offset-fs*1;
            attempted_stop_idx = IDX(i)+ttl_diode_offset+fs*1-1;
            
            %Cut epoch
            epoch = EEG_highpass_cleanline.data(:, max([1, attempted_start_idx]):min([attempted_stop_idx, length(EEG_highpass_cleanline.data)]));
            epoch = -(epoch-repmat(nanmean(epoch(:,1:fs),2),[1,size(epoch,2)]));

            % Pad with NaN if not enough data
            if attempted_start_idx < 1
                epoch = [NaN(size(epoch, 1), 1 - attempted_start_idx), epoch];
            end
            if attempted_stop_idx > length(EEG_highpass_cleanline.data)
                epoch = [epoch, NaN(size(epoch, 1), attempted_stop_idx - length(EEG_highpass_cleanline.data))];
            end
            
            epochs(:,:,i) = epoch;
            
        end
        
        %Save epochs
        cd(save_dir)
        save('face_2000ms_epochs.mat','epochs','header','-v7.3');

        %% Cut longer epoch -3 to +3
        
        %Define the time post event onset to cut epoch
        delay_identifier = 3;
        
        %Only cut longer epochs for faces
        if type == 1 
            
            epochs = [];
            
            %Loop over the IDX variable
           for i = 1:length(IDX)
                
                %Define start and stop times of epoch
                attempted_start_idx = IDX(i)+ttl_diode_offset-fs*3;
                attempted_stop_idx = IDX(i)+ttl_diode_offset+fs*delay_identifier-1;
                
                %epoch = NaN(size(EEG_highpass_cleanline, 1), fs * 8);
                
                %Cut epoch
                epoch = EEG_highpass_cleanline.data(:, max([1, attempted_start_idx]):min([attempted_stop_idx, length(EEG_highpass_cleanline.data)]));
                epoch = -(epoch-repmat(nanmean(epoch(:,1:3*fs),2),[1,size(epoch,2)]));
                
                % Pad with NaN if not enough data
                if attempted_start_idx < 1
                    epoch = [NaN(size(epoch, 1), 1 - attempted_start_idx), epoch];
                end
                
                if attempted_stop_idx > length(EEG_highpass_cleanline.data)
                    epoch = [epoch, NaN(size(epoch, 1), attempted_stop_idx - length(EEG_highpass_cleanline.data))];
                end
                
                epochs(:, 1:size(epoch, 2), i) = epoch;
                
           end
            
            %Save epochs
            cd(save_dir)
            save('face_6000ms_epochs.mat','epochs','header','-v7.3');
        
        end

    end
    
    %% Finalize TTL pulse identifiers for saving and post-processing
    
    %Find the onset times for the pulses of interest
    pulse_of_interest_sample_idx = IDX;

    %Save a matrix of the identifier (under the 'Location Accuracy'
    % column) so you can sort it later
    trial_identifiers = raw(403:end,location); %Only consider rows after calibration
  
    %Save pulse/ttl times
    cd(save_dir)
    save('ttl_times_and_identity.mat','pulse_of_interest_sample_idx', 'trial_identifiers', '-v7.3')

    %% epoch level rejections 
    CP_epochs = epochs(:,:,strcmp(trial_identifiers,'Confirmed perceived')) ;
    CnP_epochs = epochs(:,:,strcmp(trial_identifiers,'Confirmed not perceived')) ;
    
    cd('V:\RNS Study\Analysis\Analysis Code\Preprocessing\Preprocessing Pipeline Functions')
    CP_all_epochs_reref_interp = epoch_level_bad_channel_removal_RNS(CP_epochs, 'CP', save_dir, voltage_threshold, high_volt_prop, EEG_blank, prop_bad_channels);
    cd('V:\RNS Study\Analysis\Analysis Code\Preprocessing\Preprocessing Pipeline Functions')
    CnP_all_epochs_reref_interp = epoch_level_bad_channel_removal_RNS(CnP_epochs, 'CnP', save_dir, voltage_threshold, high_volt_prop, EEG_blank, prop_bad_channels);
    
    %% get channel names 
    electrode_names = string([str2mat(EEG_blank.chanlocs.labels)]); 
    preprocessed_dir = strcat('\\blumenfeld1.med.yale.internal\data7\RNS Study\Analysis\Subject Analysis\',ID,'\Scalp EEG Analysis\Preprocessed');
    if ~exist(preprocessed_dir)
        mkdir(preprocessed_dir)
    end
    cd(preprocessed_dir)
    save post_rej_epochs CP_all_epochs_reref_interp CnP_all_epochs_reref_interp electrode_names
    
end