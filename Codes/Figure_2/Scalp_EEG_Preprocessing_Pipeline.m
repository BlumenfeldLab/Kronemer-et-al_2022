%% Preprocessing Scalp EEG Data Pipeline

%This code will take raw scalp data and behavioral files and run through
%our in-house artifact rejection pipeline. 

% Written by: Sharif I. Kronemer, Mark Aksen, Julia Ding
% Date: 11/11/2019
% Modified date: 4/30/2021

clear

%% Subject ID
prompt_3 = 'Subject ID: ';
subject_ID = input(prompt_3,'s');

%Define ID value
ID = subject_ID;

%% Select run location
prompt_1 = 'Running code local or server [l, s]: ';
run_location = input(prompt_1,'s');

%% Directories

%Define directories according to where code is run (locally or on the server)
if isequal(run_location, 'l')
    
    %Add directories to Matlab path
    addpath('T:\HNCT sEEG Study\sEEG Analysis\Analysis Code') %Code folder
    addpath('T:\HNCT sEEG Study\sEEG Analysis\Analysis Code\Preprocessing\Blink')
    addpath('T:\HNCT sEEG Study\sEEG Analysis\Analysis Code\eeglab14_0_0b') %EEGLab folder
    addpath(genpath('T:\HNCT sEEG Study\sEEG Analysis\Analysis Code\eeglab14_0_0b\functions'))
    addpath('T:\HNCT sEEG Study\sEEG Analysis\Analysis Code\Extract Epochs') 
    addpath('T:\HNCT sEEG Study\sEEG Analysis\Analysis Code\Preprocessing')
    addpath('T:\HNCT sEEG Study\sEEG Analysis\Analysis Code\Preprocessing\Preprocessing Pipeline Functions')
    addpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\EEG Analysis\eeglab14_0_0b\plugins\clean_rawdata-master')

    %Define main directories
    home_dir = 'T:\HNCT sEEG Study\sEEG Analysis\Subject Analysis'; %Subject folder
    cleanline_dir = 'T:\HNCT sEEG Study\sEEG Analysis\Analysis Code\eeglab14_0_0b\plugins\tmullen-cleanline-696a7181b7d0'; %Cleanline function folder
    
    %Load blank EEGlab structure
    load('T:\HNCT sEEG Study\sEEG Analysis\Analysis Code\eeglab_prepared_blank.mat');
    
elseif isequal(run_location, 's')

    %Add directories to Matlab path
    addpath('//mnt/Data6/HNCT sEEG Study/sEEG Analysis/Analysis Code') %Code folder
    addpath('//mnt/Data6/HNCT sEEG Study/sEEG Analysis/Analysis Code/Preprocessing/Blink')
    addpath('//mnt/Data6/HNCT sEEG Study/sEEG Analysis/Analysis Code/eeglab14_0_0b') %EEGLab folder
    addpath(genpath('//mnt/Data6/HNCT sEEG Study/sEEG Analysis/Analysis Code/eeglab14_0_0b/functions'))
    addpath('//mnt/Data6/HNCT sEEG Study/sEEG Analysis/Analysis Code/Extract Epochs')
    addpath('//mnt/Data6/HNCT sEEG Study/sEEG Analysis/Analysis Code/Preprocessing')
    addpath('//mnt/Data6/HNCT sEEG Study/sEEG Analysis/Analysis Code/Preprocessing/Preprocessing Pipeline Functions')
    addpath('//mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/EEG Analysis/eeglab14_0_0b/plugins/clean_rawdata-master')
    
    %Define main directories
    home_dir = '//mnt/Data6/HNCT sEEG Study/sEEG Analysis/Subject Analysis'; %Subject folder
    cleanline_dir = '//mnt/Data6/HNCT sEEG Study/sEEG Analysis/Analysis Code/eeglab14_0_0b/plugins/tmullen-cleanline-696a7181b7d0'; %Cleanline function folder

    %Load blank EEGlab structure
    load('/mnt/Data6/HNCT sEEG Study/sEEG Analysis/Analysis Code/eeglab_prepared_blank.mat');
    
end

%% PARAMETERS

%EEG collection rate
sampling_rate = 1000;

%Final epoch duration (ms)
epoch_duration = 2000;

%Length of epoch (*2 = full epoch length)
half_epoch_duration = 2;

%Epoch removal parameters

%Bad epoch window (-200ms to +500ms pre and post event)
bad_window = 1800:2500;

%Bad sample proportion
bad_sample_prop = 0.25;

%PCA/ICA procedure

%Define the number of components to create from ICA
num_components = 10;

%% START of PREPROCESSING

%Start timer
tic

%% Define subject specific directories
if isequal(run_location, 'l') %Running locally

    %Raw Data Directory
    EEG_dir = ['T:\HNCT sEEG Study\Subject Raw Data\',ID,'\EEG Data\Perception Task'];

    %Blink Data Directory
    %Blink_data = ['T:\HNCT sEEG Study\sEEG Analysis\Subject Analysis\', ID, '\Pupillometry Analysis\Perception Task\Epochs'];

    %Photogrammetry Directory
    Photo_dir = ['S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\EEG Analysis\eeglab14_0_0b\sample_locs\GSN-HydroCel-257.sfp'];

    %Behavioral Directory 
    Beh_dir = ['T:\HNCT sEEG Study\Subject Raw Data\',ID,'\Behavioral Data\Perception Task'];

    %Save Directory
    save_dir = ['T:\HNCT sEEG Study\sEEG Analysis\Subject Analysis\',ID,'\EEG Analysis\Perception Task'];      

elseif isequal(run_location, 's') %Run on server

    %Raw Data Directory
    EEG_dir = ['//mnt/Data6/HNCT sEEG Study/Subject Raw Data/',ID,'/EEG Data/Perception Task'];

    %Blink Data Directory
    %Blink_data = ['//mnt/Data6/HNCT sEEG Study/sEEG Analysis/Subject Analysis/', ID, '/Pupillometry Analysis/Perception Task/Epochs/'];

    %Photogrammetry Directory
    Photo_dir = ['/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/EEG Analysis/eeglab14_0_0b/sample_locs/GSN-HydroCel-257.sfp'];

    %Behavioral Directory 
    Beh_dir = ['//mnt/Data6/HNCT sEEG Study/Subject Raw Data/',ID,'/Behavioral Data/Perception Task'];

    %Save Directory
    save_dir = ['//mnt/Data6/HNCT sEEG Study/sEEG Analysis/Subject Analysis/',ID,'/EEG Analysis/Perception Task'];

end

%Preprocessed Epoch Directory
preprocessed_dir = fullfile(save_dir,'Preprocessed Data');

%Power Epoch Directory
power_dir = fullfile(preprocessed_dir,'Frequency Power Epochs');

%Define the epochs directory to epochs_dir
epochs_dir = fullfile(preprocessed_dir,'Voltage Epochs');

%Behavioral times and events types
events_dir = fullfile(save_dir,'Event Times and Identifiers');

%Define the bad channel directory
bad_channel_dir = fullfile(save_dir, 'Bad Channels and Samples');

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

%% Preface - Find the number of raw files

%Note: the subject folder architecture has the run and calibration raw
%files separated except for 248VG, hence the need for unique treatment.

%Find the folders/files in the EEG_dir directory
EEG_folders = dir(EEG_dir); 

%Remove the first two rows in the EEG_files that are . and ..
EEG_folders(1:2) = []; 

%Create a cell with all the names of the raw files for EEG data
EEG_data = {};

%Set exception for 248VG that has two behavioral files but 1 EEG set file
if isequal(ID, '248VG')

    %The first raw file contains calibration
    for n = 1:size(EEG_folders,1)   
       
       %Navigate to the folder with the EEG files
       cd(fullfile(EEG_dir, EEG_folders(n).name))

       %Find all .raw files
       raw_file = dir('*.raw');

       %Create a cell of file names with .raw extension
       EEG_data{n} = raw_file.name;

    end

else

    %Loop over raw file directories
    for n = 2:size(EEG_folders,1) %Starts from 2 because the first .raw file/folder is calibration and is currently not being analyzed

       %Navigate to the folder with the EEG files
       cd(fullfile(EEG_dir, EEG_folders(n).name))

       %Create a directory of file names with .raw extension
       raw_file = dir('*.raw');

       %Create a cell of file names with .raw extension
       EEG_data{n-1} = raw_file.name;

    end

end
    
%% STEP 1: Behavioral Analysis

%Behavioral analysis to identify perception, opacity, first stimulus, and
%location information for each face stimulus used to categorize and cut 
%trials later in preprocessing.

eeg_behavioral_analysis_report_paradigm(ID, Beh_dir, save_dir, events_dir);         

%% STEP 2: EEG DIN Extraction

%Find the DIN event times from EEG file and categorize by main event
%types - faces, questions, and button presses. Correct the face DIN
%time by 42ms due to delay of stimulus display. Also, corrects for TTL
%pulse delay.

eeg_DIN_extraction_report_paradigm(ID, EEG_dir, EEG_data, EEG_folders, events_dir, run_location);

%% STEP 3: Session Level Rejection

%This function first highpass filters the session-level data, applies
%cleanline to remove 60Hz noise, uses clean_channels and clean_windows to
%find bad session-level channels and bad timepoints, interpolates over the
%bad timepoints, and average references the data while restoring the Cz
%reference channel for a total of 257 sensors.

eeg_session_level_processing_report_paradigm(ID, EEG_dir, EEG_data, EEG_folders, bad_channel_dir, preprocessed_dir, Photo_dir, events_dir, cleanline_dir, run_location);

%% STEP 3: Epoch Level Rejection

%This function will first cut the main DIN marker epochs (irrelevant faces,
%relevant faces, questions, and button presses). Bad epochs are identified
%by the number of bad samples found by clean_window in a previous
%preprocessing function. Epochs, bad epoch indices, and good sample mask
%epochs are saved.

eeg_epoch_level_processing_report_paradigm(EEG_data, events_dir, bad_channel_dir, preprocessed_dir, sampling_rate, half_epoch_duration, bad_window, bad_sample_prop)

%% STEP 5: PCA/ICA - voltage data

%Run PCA and ICA on average referenced relevant, irrelevant, button press,
%and question epochs. This process is implemented by combining all faces
%epochs in one ICA and all questions and button presses in a second ICA.

%Run face epochs
disp('Running ICA on face epochs')

%Load average reference data
cd(preprocessed_dir)
load(['Avg_ref_main_events_',num2str(half_epoch_duration*2),'s_epochs.mat'])

%Reset ICA epochs variable
ICA_epochs = [];

%Define event_type
event_type = 'faces';

%Define epochs to be included for ICA - Combine CP and CnP Trials
ICA_epochs = All_faces_epochs;

%Run ICA    
[face_epochs_components_removed_data, faces_components_removed] = eeg_pca_ica_on_voltage_report_paradigm(ID, event_type, ICA_epochs, Photo_dir, preprocessed_dir, sampling_rate, num_components, eeglab_prepared_blank);

%Save ICA component removed voltage data
cd(epochs_dir)
save face_epochs_preprocessed_components_removed.mat face_epochs_components_removed_data faces_components_removed 

%Clear variables
%clear ICA_epochs face_epochs_components_removed_data faces_components_removed 

%{
% Run question and button press epochs
disp('Running ICA on questions and button press epochs')

%Reset ICA epochs variable
ICA_epochs = [];

%Define event_type
event_type = 'questions_buttons';

%Define epochs to be included for ICA - Combine CP and CnP Trials
ICA_epochs = cat(3, All_questions_epochs, All_button_presses_epochs);

%Run ICA    
[ICA_epochs_components_removed, remove_components] = eeg_pca_ica_on_voltage_NRP(ID, event_type, ICA_epochs, Photo_dir, preprocessed_dir, sampling_rate, num_components, eeglab_prepared_blank);

%Rename variables and cut out CP and CnP epochs
questions_epochs_components_removed = ICA_epochs_components_removed(:,:,1:size(All_questions_epochs,3));
button_presses_epochs_components_removed = ICA_epochs_components_removed(:,:,size(All_questions_epochs,3)+1:size(ICA_epochs_components_removed,3));

%Check CP and CnP epochs were cut correctly
if ~isequal(size(All_questions_epochs,3), size(questions_epochs_components_removed ,3)) || ~isequal(size(All_button_presses_epochs,3), size(button_presses_epochs_components_removed,3))

    %Warning about mismatch in the number of trials
    warning('The number of CP/CnP trials post-ICA does not match!')

end

%Rename ICA components removed
questions_button_presses_components_removed = remove_components;

%Save ICA component removed voltage data
cd(epochs_dir)
save questions_buttons_epochs_preprocessed_components_removed.mat questions_epochs_components_removed button_presses_epochs_components_removed questions_button_presses_components_removed

%Clear variables
clear ICA_epochs ICA_epochs_components_removed All_questions_epochs All_button_presses_epochs ...
    questions_epochs_components_removed button_presses_epochs_components_removed questions_button_presses_components_removed
%}
 
%% STEP 4: Cut Epochs

%This function takes the ICA component removed epochs and cuts specfic
%epoch types using the behavioral identifier information, including CP,
%CnP, and irrelevant stimuli by opacity type. This epochs are saved in the
%epochs directory. Bad epochs identified earlier for each main event type
%are excluded. 

cut_eeg_specific_epochs_report_paradigm(EEG_data, events_dir, epochs_dir, bad_channel_dir, half_epoch_duration)
    
%% STEP 7: Extract Frequency Content of Signal - Spectrogram

%Extract the frequenc content from those signals using the spectrogram methodology. Power is
%calculated using both a log and z-score approach and both are saved for
%each epoch type. Note: for computational efficiency purposes, only the
%first 50 button presses are considered. A specific invesitgation of
%button press epochs will require rerunning this frequency extraction
%with all the button press trials included.

%Data types
data_types = {'CP_1s','CP_15s','CnP_1s','CnP_15s'};

%Load bad trial index
cd(events_dir)
load('Bad_trials_index.mat')

%Loop over data types
for type = 1:size(data_types,2)

    %Define epoch type
    epoch_type = char(data_types(type));

    disp(['**Running spectrogram ', epoch_type, ' epochs'])

    %Define epochs to be included for ICA  
    cd(epochs_dir)
    load([epoch_type,'_4s_epochs.mat'])

    %Remove bad epochs using bad_epochs_idx if button press epochs     
    if isequal(epoch_type, 'button_presses')

        disp('Removing bad button press trials')

        %Remove bad epochs using bad_epochs_idx
        eval([epoch_type,'_epochs(:,:,',epoch_type,'_bad_epochs_idx) = [];'])

    end

    %Rename epochs
    epochs = eval([epoch_type,'_epochs']);

    %Cut down the number of button press epochs for computational
    %efficiency
    if isequal(epoch_type, 'button_presses')

        %Threshold of 50 button press trials
        if size(epochs,3) > 100

            %Select the first 50 button presses
            epochs = epochs(:,:,1:100);

        end

    end

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
    [electrode_all_freq_log_power_vector, electrode_all_freq_zscore_power_vector] = eeg_spectrogram_on_trials_report_paradigm(epochs, epoch_type, sampling_rate, epoch_duration);

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

eeg_preprocessing_summary_txt_file_report_paradigm(ID, EEG_data, bad_channel_dir, preprocessed_dir)

%End timer
toc
