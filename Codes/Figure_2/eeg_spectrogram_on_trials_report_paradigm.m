%% Scalp EEG Time-frequency analysis for all epoch types

% Extract frequency power from epcoh voltage using spectrogram approach. 

%Adapted/Written by: Mark Aksen and Sharif I. Kronemer
%Date: 01/30/2019
%Modified: 1/19/2021

function [electrode_all_freq_log_power_vector, electrode_all_freq_zscore_power_vector] = eeg_spectrogram_on_trials_report_paradigm(epochs, epoch_type, sampling_rate, epoch_duration)

    %Cut frequency range due to storage and computational efficiency
    freq_cut = [1:200]; %1 to 200Hz

    %Cut epochs duration (cut extra samples)
    %epochs = epochs(:,abs((epoch_duration/2)-2001):(epoch_duration/2)+2001,:);

    %Define the size of the window for each bin of power analysis
    window_size = floor(sampling_rate/8); %250ms

    %Define the number of samples you want to slide the window by
    sliding_window = floor(window_size/16); %31.25ms

    %Calculate the total number of bins that will result from your window
    epoch_size = size(epochs,2); %Find the number of samples in the epoch (e.g., 1000ms or 5000ms at 1000Hz)
    nonoverlap_bins_num = floor(epoch_size/sliding_window); %Divide epoch samples by sliding window
    bins_per_window = floor(window_size/sliding_window); %Total number of bins with a single window
    total_bins_num = nonoverlap_bins_num - bins_per_window + 1; %Calculate the total number of bins in the epoch

    %Find the number of electrodes (e.g., 256)
    nElectrode = size(epochs,1);

    %Find the number of trial
    num_trials = size(epochs,3);

    %% Time-frequency Analysis
    
    %Initialize Variable Trials(electrodes x frequency x num bin, trials)
    electrode_all_freq_log_power_vector = [];
    electrode_all_freq_zscore_power_vector = [];

    %Loop over electrodes
    for electrode = 1:nElectrode

        %Check for NaN - skip channel if NaN
        if isnan(epochs(electrode,:,1))
            
            disp(['Skipping NaN channel - ', num2str(electrode)])
            continue
            
        end   

        %Loop over trials
        for trial_num = 1:num_trials

            %Select voltage data from a single trial and electrode
            trial = squeeze(epochs(electrode, :, trial_num));

            %Run Spectrogram on single trial and electrode
            [trial_power, ~, ~] = spectrogram(trial, window_size, window_size-sliding_window, floor(sampling_rate), sampling_rate);

            %Cut down the frequency range
            trial_power = trial_power(freq_cut,:);
            
            %Calculate power - squared, absolute value
            trial_power = abs(trial_power).^2;

            %Take logarithm of power (channel x time)
            log_trial_power = log(trial_power);
            
            %Store all frequencies (electrode x frequency x time x trials)
            electrode_all_freq_log_power_vector(electrode,:,:,trial_num) = log_trial_power; %Log power
            electrode_all_freq_zscore_power_vector(electrode,:,:,trial_num) = trial_power; %Non-log power

        end

    end

    %% Baseline correction (subtract mean power across each electrode individually)
    
    disp('Running log power')
    
    %Loop over frequencies
    for freq = 1:size(electrode_all_freq_log_power_vector,2)

        %Loop over electrodes
        for electrode = 1:nElectrode

            %Loop over trials
            for trial_num = 1:num_trials

                %Baseline power - baseline equals mean power in prestimulus period
                electrode_all_freq_log_power_vector(electrode,freq,:,trial_num) = electrode_all_freq_log_power_vector(electrode,freq,:,trial_num)...
                    -nanmean(electrode_all_freq_log_power_vector(electrode,freq,1:floor(total_bins_num/2),trial_num));

            end

        end

    end

    %% Zscore Non-log Power

    disp('Running z-score power')
    
    %Loop over frequencies
    for freq = 1:size(electrode_all_freq_zscore_power_vector,2)

        %Loop over electrodes
        for electrode = 1:nElectrode

            %Find the mean baseline and standard deviation for all
            %trials of that reference range for that electrode -
            %dimensions channel x frequency x time x trials
            all_trial_mean = nanmean(nanmean(squeeze(electrode_all_freq_zscore_power_vector(electrode,freq,1:floor(total_bins_num/2),:)),2),1);
            all_trial_SD = std(nanmean(squeeze(electrode_all_freq_zscore_power_vector(electrode,freq,1:floor(total_bins_num/2),:)),2),[],1);

            %Loop over trials
            for trial_num = 1:num_trials

                %Define freq, trial, and electrode to z-score
                data = squeeze(electrode_all_freq_zscore_power_vector(electrode,freq,:,trial_num));

                %Baseline power - baseline equals mean power in prestimulus period
                electrode_all_freq_zscore_power_vector(electrode,freq,:,trial_num) = (data-all_trial_mean)/all_trial_SD;

            end

        end

    end
    
    %% Average over trials for button presses to save space
    
    %Note: This can be used when running all button press events that will
    %generate a large variable (10GB+) if all the trails are saved instead
    %of average over the trials before saving. 
    
%     %If the event is button presses
%     if isequal(epoch_type, 'button_presses')
%         
%         %Final matrix is 3D - channel x frequency x time
%         electrode_all_freq_log_power_vector = mean(electrode_all_freq_log_power_vector, 4);
%         electrode_all_freq_zscore_power_vector = mean(electrode_all_freq_zscore_power_vector, 4);
%     
%     end
    
end
