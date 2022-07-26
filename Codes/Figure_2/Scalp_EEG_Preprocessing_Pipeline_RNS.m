%% Preprocessing Scalp EEG Data Pipeline

%This code will take raw scalp data and behavioral files and run through
%our in-house artifact rejection pipeline. 

%Step 1: Behavioral analysis
%Step 2: EEG TTL time extraction - Computed with seperate script
%Step 3: Session level rejection and processing
%Step 4: Epoch level rejection
%Step 5: PCA/ICA
%Step 6: Cut specfic epochs
%Step 7: Frequency extraction
%Step 8: Create Preprocessign Summary Report

% Written by: Sharif I. Kronemer, Mark Aksen, & Julia Ding
% Date: 11/11/2019
% Modified date: 6/4/2021

clear

%% Prompts

% Subject ID
prompt_1 = 'Subject ID: ';
subject_ID = input(prompt_1,'s');

%Define ID value
ID = subject_ID;

% Select run location
prompt_2 = 'Running code local or server [l, s]: ';
run_location = input(prompt_2,'s');

% Reference type
prompt_3 = 'Reference type [average, mastoid]: ';
reference_type = input(prompt_3,'s');

%Define reference names
if isequal(reference_type,'average')
    
    reference_name = 'Average Reference';
    
elseif isequal(reference_type,'mastoid')
    
    reference_name = 'Mastoid Reference';
    
end

%% Directories

%Define directories according to where code is run (locally or on the server)
if isequal(run_location, 'l')
    
    %Add directories to Matlab path
    addpath('\\server2.med.yale.internal\Data8\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code') %Code folder
    addpath('\\server2.med.yale.internal\Data8\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\EEG Analysis\eeglab14_0_0b') %EEGLab folder
    addpath(genpath('\\server2.med.yale.internal\Data6\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\EEG Analysis\eeglab14_0_0b\functions'))
    addpath('\\server2.med.yale.internal\Data8\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\EEG Analysis\eeglab14_0_0b\plugins\clean_rawdata-master')
    addpath('\\server2.med.yale.internal\Data7\RNS Study\Analysis\Analysis Code\Epoch extraction functions')
    addpath('\\server2.med.yale.internal\Data7\RNS Study\Analysis\Analysis Code\Preprocessing\Scalp EEG Preprocessing\Preprocessing Pipeline Functions')
    
    %Define main directories
    home_dir = '\\server2.med.yale.internal\Data7\RNS Study\Analysis\Subject Analysis'; %Subject folder
    cleanline_dir = '\\server2.med.yale.internal\Data8\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\EEG Analysis\eeglab14_0_0b\plugins\tmullen-cleanline-696a7181b7d0'; %Cleanline function folder
    
    %Load blank EEGlab structure
    load('\\server2.med.yale.internal\Data7\RNS Study\Analysis\Analysis Code\Preprocessing\Scalp EEG Preprocessing\RNS_eeglab_prepared_blank.mat');

elseif isequal(run_location, 's')

    %Add directories to Matlab path
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code') %Code folder
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/EEG Analysis/eeglab14_0_0b') %EEGLab folder
    addpath(genpath('/mnt/Data6/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/EEG Analysis/eeglab14_0_0b/functions'))
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/EEG Analysis/eeglab14_0_0b/plugins/clean_rawdata-master')
    addpath('/mnt/Data7/RNS Study/Analysis/Analysis Code/Epoch extraction functions')
    addpath('/mnt/Data7/RNS Study/Analysis/Analysis Code/Preprocessing/Scalp EEG Preprocessing/Preprocessing Pipeline Functions')
    
    %Define main directories
    home_dir = '/mnt/Data7/RNS Study/Analysis/Subject Analysis'; %Subject folder
    cleanline_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/EEG Analysis/eeglab14_0_0b/plugins/tmullen-cleanline-696a7181b7d0'; %Cleanline function folder
    
    %Load blank EEGlab structure
    load('/mnt/Data7/RNS Study/Analysis/Analysis Code/Preprocessing/Scalp EEG Preprocessing/RNS_eeglab_prepared_blank.mat');
    
end

%Store electrode names
electrode_names = [];

%Loop over channels
for channel = 1:20
    
    %Store electrode name
    electrode_names{channel,1} = char(RNS_eeglab_prepared_blank.chanlocs(channel).labels);
    
end

%% PREPROCESSING PARAMETERS

%EEG collection rate
sampling_rate = 256;

%Final epoch duration (ms)
epoch_duration = 2000;

%Length of epoch (*2 = full epoch length)
half_epoch_duration = 2;

%Epoch removal parameters

%Bad epoch window (-200ms to +500ms pre and post event) Note: 256 sampling
%rate makes the face time 512ms
bad_window = 461:640; %Could be written as round(1.8*sampling_rate):round(2.5*sampling_rate)

%Bad sample proportion
bad_sample_prop = 0.25;

%PCA/ICA procedure

%Define the number of components to create from ICA
num_components = 10;

%% START of PREPROCESSING

%Start timer
tic

%% Define subject specific directories

if isequal(run_location, 'l') %Run locally

    %Raw Data Directory
    EEG_dir = ['\\server2.med.yale.internal\Data7\RNS Study\Subject Raw Data\',ID,'\Scalp EEG Data\Calibration and All Runs'];

    %Photogrammetry Directory
    Photo_dir = ['\\server2.med.yale.internal\Data8\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\EEG Analysis\eeglab14_0_0b\sample_locs\Standard-10-20-Cap81.locs'];

    %Behavioral Directory 
    Beh_dir = ['\\server2.med.yale.internal\Data7\RNS Study\Subject Raw Data\',ID,'\Behavioral Data\Calibration and All Runs'];

    %Save Directory
    save_dir = ['\\server2.med.yale.internal\Data7\RNS Study\Analysis\Subject Analysis\',ID,'\Scalp EEG Analysis'];

elseif isequal(run_location, 's') %Run on server

    %Raw Data Directory
    EEG_dir = ['/mnt/Data7/RNS Study/Subject Raw Data/',ID,'/Scalp EEG Data/Calibration and All Runs'];

    %Photogrammetry Directory
    Photo_dir = ['/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/EEG Analysis/eeglab14_0_0b/sample_locs/Standard-10-20-Cap81.locs'];

    %Behavioral Directory 
    Beh_dir = ['/mnt/Data7/RNS Study/Subject Raw Data/',ID,'/Behavioral Data/Calibration and All Runs'];

    %Save Directory
    save_dir = ['/mnt/Data7/RNS Study/Analysis/Subject Analysis/',ID,'/Scalp EEG Analysis'];

end

%Preprocessed Epoch Directory
preprocessed_dir = fullfile(save_dir,'Preprocessed Data',reference_name);

%Power Epoch Directory
power_dir = fullfile(preprocessed_dir,'Frequency Power Epochs');

%Define the epochs directory to epochs_dir
epochs_dir = fullfile(preprocessed_dir,'Voltage Epochs');

%Behavioral times and event types
events_dir = fullfile(save_dir,'Event Times and Identifiers');

%Define the bad channel directory
bad_channel_dir = fullfile(save_dir,'Bad Channels and Samples',reference_name);

%Define blink directory 
blink_index_dir = fullfile(save_dir,'Bad Channels and Samples');

%Check for the existences of the epochs directory
if ~exist(epochs_dir)
    
    %Create folder
    mkdir(epochs_dir); 

end    

%Check for the existence of the events directory
if ~exist(events_dir)

    %Create folder
    mkdir(events_dir)

end

%Check for the existence of the bad channel directory
if ~exist(bad_channel_dir)

    %Create folder
    mkdir(bad_channel_dir)

end

%Check for the existence of the preprocessed directory
if ~exist(preprocessed_dir)

    %Create folder
    mkdir(preprocessed_dir);

end

%Check for the existence of the power directory
if ~exist(power_dir)

    %Create folder
    mkdir(power_dir);

end

%Open EEGLab without GUI
eeglab('nogui');

%% STEP 1: Behavioral Analysis

%Behavioral analysis to identify perception and hemifield. 

eeg_behavioral_analysis_RNS(ID, Beh_dir, events_dir);

%% STEP 2: EEG TTL Pulse Extraction

%TTL pulse extraction implemented in a separate preprocessing scripts
%(e.g., preprocessing_536JD.m) - These scripts are specially designed for
%each subject to manage some of the artifacts and errors in TTL pulse
%identification


%% STEP 3: Session Level Rejection

%This function first highpass filters the session-level data, applies
%cleanline to remove 60Hz noise, uses clean_channels and clean_windows to
%find bad session-level channels and bad timepoints, and average references the data.

eeg_session_level_processing_RNS(ID, events_dir, EEG_dir, bad_channel_dir, preprocessed_dir, Photo_dir, cleanline_dir, RNS_eeglab_prepared_blank, reference_type);

%% STEP 4: Epoch Level Rejection

%This function will first cut the main face epochs by TTL pulse time. Bad epochs are identified
%by the number of bad samples found by clean_window in a previous
%preprocessing function. Epochs, bad epoch indices, and good sample mask
%epochs are saved.

eeg_epoch_level_processing_RNS(ID, EEG_dir, events_dir, bad_channel_dir, preprocessed_dir, sampling_rate, half_epoch_duration, bad_window, bad_sample_prop);

%% STEP 5: PCA/ICA - voltage data

%Run PCA and ICA on average referenced face epochs. 

%Load average reference data
cd(preprocessed_dir)
load(['Main_events_',num2str(half_epoch_duration*2),'s_epochs.mat'])

%Run face epochs
disp('Running ICA on face epochs')

%Reset ICA epochs variable
ICA_epochs = [];

%Define event_type
event_type = 'faces';

%Define epochs to be included for ICA
ICA_epochs = All_faces_epochs;

%Run ICA    
[face_epochs_components_removed, face_components_removed] = eeg_pca_ica_on_voltage_RNS(event_type, ICA_epochs, Photo_dir, preprocessed_dir, sampling_rate, num_components, RNS_eeglab_prepared_blank);

%Save ICA component removed voltage data
cd(epochs_dir)
save face_epochs_preprocessed_components_removed.mat face_epochs_components_removed face_components_removed

%Clear variables
clear ICA_epochs ICA_epochs_components_removed All_faces_epochs ...
    faces_epochs_components_removed faces_components_removed

%% STEP 5: TTL Pulse Artifact Rejection

%Average referencing note: Common avg referencing eliminates the TTL pulse
%artifact because it is found across all channels.

%Note: All subjects to varying degrees show an signature TTL pulse artifact
%around the time of stimulus onset

%ICA for artifact parameters 
num_components = 5; %471MH (remove #1), 489RD(remove #2), 490KB (remove #2), 491GS (remove #1), 535BP(removed #1), 536JD did not have TTL pulses
event_type = 'TTL';
TTL_window = 460:580; %512 is stimulus onset, 256 samples rate, epoch length 1025 data points (4 seconds) 

%Load data
cd(epochs_dir)
load('face_epochs_preprocessed_components_removed.mat')

%Cut data to TTL pulse artifact time
ICA_epochs = face_epochs_components_removed(:,TTL_window,:);

%Run ICA    
[ICA_epochs_components_removed, TTL_components_removed] = eeg_pca_ica_on_voltage_RNS(event_type, ICA_epochs, Photo_dir, preprocessed_dir, sampling_rate, num_components, RNS_eeglab_prepared_blank);

%Replace TTL time in original data with TTL ICA removed data
face_epochs_components_removed(:,TTL_window,:) = ICA_epochs_components_removed;

%Save data and components
cd(epochs_dir)
save face_epochs_preprocessed_ttl_components_removed.mat face_epochs_components_removed TTL_components_removed

%% STEP 6: Cut Specifc Epochs - CP, CnP, Irrelevant Opaque, etc.

%This function takes the ICA component removed epochs and cuts specfic
%epoch types using the behavioral identifier information. This epochs are saved in the
%epochs directory. Bad epochs identified earlier for each main event type
%are excluded. 

cut_eeg_specific_epochs_RNS(events_dir, epochs_dir, bad_channel_dir, half_epoch_duration, blink_index_dir)

%% STEP 7: Extract Frequency Content of Signal - Spectrogram

%Extract the frequenc content from those signals using the spectrogram methodology. Power is
%calculated using both a log and z-score approach and both are saved for
%each epoch type. Output is a 4D matrix with the dimensions: channels x
%frequency x time bins x trials. 

%Load CP and CnP epochs 
cd(epochs_dir)
load('CP_4s_epochs.mat')
load('CnP_4s_epochs.mat')

%Data types
data_types = {'CP_good', 'CP_all', 'CnP_good', 'CnP_all'};

%Loop over data types
for type = 1:size(data_types,2)

    %Define epoch type
    epoch_type = char(data_types(type));

    disp(['Running spectrogram ', epoch_type, ' epochs'])

    %Rename epochs
    epochs = eval([epoch_type,'_epochs']);

    %Check if there are no trials of this event type
    if isempty(epochs)

       %Create an empty vector for the missing trial type to be saved later
       eval([epoch_type,'_spec_log_power_epochs = [];'])
       eval([epoch_type,'_spec_z_power_epochs = [];'])

       %Save Epochs
       cd(power_dir)
       save([epoch_type,'_spec_power_epochs.mat'], [epoch_type,'_spec*'], '-v7.3')

       %Skip to next event type
       continue 

    end

    %Run log spectrogram - 4D output: channel x frequencies, time, trials
    [electrode_all_freq_log_power_vector, electrode_all_freq_zscore_power_vector] = eeg_spectrogram_on_trials_RNS(epochs, sampling_rate, epoch_duration);

    %Rename variable with data type
    eval([epoch_type,'_spec_log_power_epochs = electrode_all_freq_log_power_vector;'])
    eval([epoch_type,'_spec_z_power_epochs = electrode_all_freq_zscore_power_vector;'])

    %Save Epochs
    cd(power_dir)
    save([epoch_type,'_spec_power_epochs.mat'], [epoch_type,'_spec*'], '-v7.3')

    %Clear Variables
    eval(['clear ',epoch_type,'_spec*'])
    clear electrode_all_freq_log_power_vector electrode_all_freq_zscore_power_vector

end

%% STEP 8: Create a Summary Report of Preprocessing

%This code will create a text files that summarizes the preprocessing
%rejection stats - session level and epoch level rejections

eeg_preprocessing_summary_txt_file_RNS(ID, bad_channel_dir, preprocessed_dir, blink_index_dir)

%End timer
toc
