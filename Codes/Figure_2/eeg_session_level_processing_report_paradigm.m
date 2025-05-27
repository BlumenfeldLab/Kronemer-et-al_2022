%% Process session level EEG data

%***Inputs***
%Subject ID
%Behavioral data directory
%EEG data directory
%Save directory
%Photogrammetry directory
%Behavioral results directory
%Cleanline directory
%Run location information

%***Outputs***
%No outputs (data is saved within function)

%Written by: Sharif I. Kronemer and Mark Aksen
%Date: 7/29/2019
%Edited: 4/30/2021

%% Loop through all the raw files available for testing session
 
function eeg_session_level_processing_report_paradigm(ID, EEG_dir, EEG_data, EEG_folders, bad_channel_dir, preprocessed_dir, Photo_dir, events_dir, cleanline_dir, run_location)
   
    %Add plugins folder to path (Note: there are conflicting functions in
    %EEGLab and Matlab so be careful not to add all EEGLab directories to
    %the path)
    if isequal(run_location, 'l') %Local
        
        addpath(genpath('T:\HNCT sEEG Study\sEEG Analysis\Analysis Code\eeglab14_0_0b\functions'))
        
    else isequal(run_location, 's') %Server
        
        addpath(genpath('mnt/Data6/HNCT sEEG Study/sEEG Analysis/Analysis Code/eeglab14_0_0b/functions'))
        
    end

    %Loop through the raw_files that equals the number of available behavioral run csv files
    for raw_files = 1:size(EEG_data,2)

        disp(['Running EEG analysis for session ', num2str(raw_files+1)])
        
        %Load variable of interest
        cd(events_dir)
        load(['Raw_file_',num2str(raw_files+1),'_trial_identifiers.mat'])

        %% EEG Data .raw Files
        
        disp('Load EEG data')

        %Use the pop_readegi function to open a raw file
        
        if strcmp(ID,'248VG')
            %EEG = pop_readegi('T:\HNCT sEEG Study\Subject Raw Data\248VG\EEG Data\Perception Task\Calibration Run 1 and 2\248VG_Visual_Task_11_2_17_1.raw');
            EEG = pop_readegi(fullfile(EEG_dir, EEG_folders(raw_files).name, EEG_data{raw_files}));

        else 
            EEG = pop_readegi(fullfile(EEG_dir, EEG_folders(raw_files+1).name, EEG_data{raw_files}));
            
        end
        
        %Error of data struct is empty
        if isempty(EEG)
            
            error('EEG structure is empty')
            
        end 
        
        %% High pass filtering of the EEG data (passband frequency = 0.1 Hz)
        
        %Note: clean_drift is a EEGLab function that utilizes filtfilt
        %procedure to highpass filter the session-level EEG data.
        
        %Transition band 
        highpass_band = [0.25 0.75];
        
        disp('1Hz Highpass Data')
        highpassEEG = clean_drifts(EEG,highpass_band);
        
        %% Applying cleanline (removing 60Hz noise)
        
        disp('Applying cleanline to remove 60 Hz noise')
        
        %Add plugins folder of EEGlab to path
        addpath(genpath(cleanline_dir))
        
        %Find the channel locations - saving this sensor location information
        %because cleanline removes this information
        chanlocs = highpassEEG.chanlocs;
        
        %Run cleanline
        cleanlineEEG = pop_cleanline(highpassEEG, 'bandwidth', 2,'chanlist', [1:highpassEEG.nbchan], 'computepower', 0, 'linefreqs', [60 120],...
        'normSpectrum', 0, 'p', 0.01, 'pad', 2, 'plotfigures', 0, 'scanforlines', 1, 'sigtype', 'Channels', 'tau', 100,...
        'verb', 1, 'winsize', 4, 'winstep', 2);
        
        %Define the channel locations in the EEG structure
        cleanlineEEG.chanlocs = chanlocs;
        
        %Remove plugins folder of EEGlab from path because it conflicts
        %with other functions
        rmpath(genpath(cleanline_dir))
        
        %% Session-level Channel Rejection
        
        %Parameters
        chancorr_crit = 0.8;
        channel_crit_maxbad_time = 0.5;
        line_crit = 4;
        
        disp('Session-level channel rejection')
        [~,bad_channels] = clean_channels(cleanlineEEG,chancorr_crit,line_crit,[],channel_crit_maxbad_time); 
        
        %Define electrode number
        bad_channels = find(bad_channels);
        
        %Rename variable
        eval(['Session_', num2str(raw_files), '_bad_channels = bad_channels;']);

        %% Session-level Sample Rejection
        
        %Parameters
        window_crit_tolerances = [-Inf,7]; %Min and max SD %NOTE: the low end is set to negative 
        %infinity so that rejections aren't made because their too quiet; While
        %the default value in clean_windows is -3.5, 5; default in
        %clean_artifact and pop_clean_rawdata is -Inf, 7;
        window_crit = 0.25; %Maximum proportion of bad channels
        
        disp('Reject time periods')
        
        %Bad samples is a mask of good (1) and bad (0) times periods
        [cleanwindowEEG, bad_samples] = clean_windows(cleanlineEEG,window_crit,window_crit_tolerances);
        
        %Visualize rejected EEG data - new data vs old data
        vis_artifacts(cleanwindowEEG, cleanlineEEG);
        
        %Review data (Need to type into command line)
%         prompt_review = 'Review data. Continue or break? [c/b]: ';
%         review = input(prompt_review,'s');
    
        review = 'c' ;
        %Break
        if isequal(review, 'b')
            
            error('Prompted by user to break from sample rejection')
            
        end  
        
        %Close figure
        close all
        
        %Check bad_samples has the correct number of timepoints 
        if ~isequal(length(bad_samples),size(cleanlineEEG.data,2))
        
            error('Problem with bad_sample size!')
        
        end
    
        %Rename variable
        eval(['Session_', num2str(raw_files), '_bad_samples = bad_samples;']);

        %% Interpolate over bad channels

        disp('Spherical interpolation')
        
        %Set channel locations
        cleanlineEEG.chanlocs = readlocs(Photo_dir); 
        cleanlineEEG.urchanlocs = readlocs(Photo_dir);
        
        %Remove Cz electrode (Note: if you don't remove it here pop_interp
        %will find a mismatch in the number of channel locations and the
        %data structure. Cz is added back in average referencing to correct
        %this mistmatch)
        cleanlineEEG.chanlocs(257) = []; 
        cleanlineEEG.urchanlocs(257) = [];

        %Interpolation - Input: dataset, channels to interpolate, and method)
        interpEEG = pop_interp(cleanlineEEG, bad_channels, 'spherical');

        %% Average Reference

        disp('Average reference signal')

        %Set channel locations
        interpEEG.chanlocs = readlocs(Photo_dir); 
        interpEEG.urchanlocs = readlocs(Photo_dir);

        %Average reference - Interpolated EEG data, average reference [], reference
        %electrode info (CZ), and all electrode location info
        EEG_avgref = reref(interpEEG.data,[],'refloc',interpEEG.chanlocs(257),'elocs',interpEEG.chanlocs);
        
        full_257_chanlocs = interpEEG.chanlocs ; 
        full_257_urchanlocs = interpEEG.urchanlocs; 
        save(fullfile(preprocessed_dir,'full_257_locs.mat'),'full_257*')
        
        %Rename data 
        eval(['Session_', num2str(raw_files), '_EEG_avgref = EEG_avgref;']);

    end
    
    %Save preprocessed session average referenced data
    cd(preprocessed_dir)
    save Session_preprocessed_data.mat Session*EEG_avgref -v7.3
    
    %Save index of bad channels and samples
    cd(bad_channel_dir)
    save Session_bad_channels_samples.mat Session*bad_channels Session*bad_samples -v7.3
 
end