%% Process Session-Level EEG data

%This function first highpass filters the session-level data, applies
%cleanline to remove 60Hz noise, uses clean_channels and clean_windows to
%find bad session-level channels and bad timepoints, interpolates over the
%bad timepoints, and average references the data while restoring the Cz
%reference channel for a total of 257 sensors.

%Written by: Sharif I. Kronemer and Mark Aksen
%Date: 7/29/2019
%Modified: 2/3/2021
 
function eeg_session_level_processing_NRP(ID, Beh_dir, events_dir, EEG_dir, bad_channel_dir, preprocessed_dir, Photo_dir, cleanline_dir)

    %Open the .raw file you intend to analyze (session #)
    clear EEG
    
    %Create structure with the names of the raw data files across runs
    cd(Beh_dir)
    Beh_data = dir('*.csv'); 
    
    %Find the folders/files in the EEG_dir directory
    cd(EEG_dir)
    EEG_data = dir([ID,'*.raw']); 
    
    %Check if number of behavioral files equals number of raw files
    if ~isequal(size(Beh_data,1),size(EEG_data,1))
        
        warning('***The number of behavioral files does not equal the number of EEG files!***')
        
    end

    %Loop through the raw_files that equals the number of available behavioral run csv files
    for raw_files = 1:size(EEG_data,1)

        disp(['Running EEG analysis - Session ', num2str(raw_files)])
        
        %Load variable of interest
        cd(events_dir)
        load(['Raw_file_',num2str(raw_files),'_trial_identifiers.mat'])

        %% EEG Data .raw Files
        
        disp('Load EEG data')

        %Use the pop_readegi function to open a raw file
        EEG = pop_readegi(fullfile(EEG_dir, EEG_data(raw_files).name));
        
        %If EEG variable is empy throw error
        if isempty(EEG) 
            
            error('EEG structure is empty') 
            
        end 
          
        %% Highpass filter 1Hz 
        
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
        
        %Select which components to remove (Need to type into command line)
        prompt_review = 'Review data. Continue or break? [c/b]: ';
        review = input(prompt_review,'s');
        
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
        
        %{
        
        %NOTE: Old artifact rejection code that also creates figures for
        the distribution of channels with high standard deviation. These
        processes have been replaced with clean_channels and clean_windows.
        
        %% Artifactual channel identification - High SD Electrodes 

        disp('Apply high SD detection')

        %Calculate and find high index channels
        [high_SD_index] = Find_high_SD_electrodes_NRP(EEG.data, ...
            EEG.event(1).latency, EEG.event(end).latency, SD_threshold, ID, raw_files, bad_channel_dir);

        %% Artifactual channel identification - High high-frequency power

        disp('Apply high frequency power detection')
        
        %Initialize parameters
        pass_frequency = 30; %Test frequency
        
        noise_threshold = 2; %Z-score high frequency power threshold
        noisy_windows_prop = 0.25; %Threshold proportion of noisy bins
        
        window_size = 4; %In seconds
        window_overlap = 2; %In seconds

        %Create filter parameters (30Hz + highpass)
        hpFilt = designfilt('highpassiir','FilterOrder',8, ...
        'PassbandFrequency',pass_frequency,'PassbandRipple',0.2, 'SampleRate',sampling_rate);

        %Filter each channel
        for c = EEG.nbchan:-1:1

            %Implement filter
            filt_data(:,c) = filtfilt(hpFilt,double(EEG.data(c,:)')); %Edited Julia

        end

        %Bin the filtered data
        nWindows = floor(EEG.pnts/(window_size*sampling_rate));
        noisyChannelWindows = zeros(EEG.nbchan, nWindows);

        %Looping over bins
        for i = 1:nWindows

            %Find the earliest and latest samples in a bin
            window_start = ((window_size-window_overlap)*sampling_rate)*(i-1)+1; 
            window_end = min(window_start + window_size*sampling_rate - 1,EEG.pnts);

            %Find the power of the filtered signal and average across
            %time for each channel 
            noisiness = mean(filt_data(window_start:window_end,:).^2,1);

            %Zscore across channels
            znoise = zscore(noisiness);

            %Create index of bins that have a zscore high frequency
            %power greater than noise_threshold
            noisy_chans_index = znoise > noise_threshold;
            noisyChannelWindows(noisy_chans_index,i) = 1;

        end

        %Clear variable
        clear filt_data

        %Find the proportion of bad bins per channel (average across
        %bins of the 0/1 index)
        noisinessMean = mean(noisyChannelWindows,2);

        %Identify channels with more than noisy_windows_prop proportion
        %of bad bins
        noisy_chans_index = find(noisinessMean > noisy_windows_prop);
        
        %Plot histogram of SD 
        histo = figure;
        hold on

        title([ID, ' Session ',num2str(raw_files),' - +30Hz Power Across All Electrodes'])
        ylabel('Number of Electrodes')
        xlabel('Proportion of Bins w/ >2 Z-score +30Hz Power')
        xlim([0 1])

        histogram(noisinessMean)
        plot([noisy_windows_prop,noisy_windows_prop],[0,100], 'r');

        %Save figure
        cd(bad_channel_dir)
        savefig([ID,'_session_',num2str(raw_files),'_30Hz_power_histogram.fig'])
        saveas(histo,[ID,'_session_',num2str(raw_files),'_3Hz_power_histogram.jpg'])

        %% Aggregate bad channels from the previous artifact procedures

        %Create variable of bad channels
        eval(['Session_', num2str(raw_files), '_bad_channels_highSD = high_SD_index;']);
        eval(['Session_', num2str(raw_files), '_bad_channels_highfreq30 = noisy_chans_index;']);   
        eval(['Session_', num2str(raw_files), '_bad_channels_all = union(high_SD_index, noisy_chans_index);']);
       %} 
        
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