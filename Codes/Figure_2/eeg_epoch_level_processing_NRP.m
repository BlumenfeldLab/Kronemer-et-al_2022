%% Epoch level bad channel removal

%This function will first cut the main DIN marker epochs (irrelevant faces,
%relevant faces, questions, and button presses). Bad epochs are identified
%by the number of bad samples found by clean_window in a previous
%preprocessing function. Epochs, bad epoch indices, and good sample mask
%epochs are saved.

%Written by: Sharif I. Kronemer & Mark Aksen
%Date: 7/29/2019
%Modified date: 5/30/2022

function eeg_epoch_level_processing_NRP(ID, EEG_dir, events_dir, bad_channel_dir, preprocessed_dir, sampling_rate, relevant_location, half_epoch_duration, bad_window, bad_sample_prop)

    %Find the folders/files in the EEG_dir directory
    cd(EEG_dir)
    EEG_data = dir([ID,'*.raw']); 

    %Load preprocessed session level data
    cd(preprocessed_dir)
    load('Session_preprocessed_data.mat')

    %Load session level bad channels/samples
    cd(bad_channel_dir)
    load('Session_bad_channels_samples.mat')

    %Main epoch types
    epoch_types = {'relevant_faces_post_jitter','irrelevant_faces_post_jitter'}; %{'relevant_faces', 'irrelevant_faces', 'button_presses', 'questions'};

    %% Cut Epochs for Face, Question, and Button Press Times
    % These lines of code will cut the epochs for question and button
    % presss events. In addition, face events will be cut, however the face
    % events are not first categorized (i.e., CP, CnP, etc.). Therefore
    % these face epochs represent all TTL pulses corresponding with a face
    % event, including blank trials. 

    %Loop over raw files
    for raw_files = 1:size(EEG_data,1)

        disp(['**Session ', num2str(raw_files), ' - Cutting epochs'])
        
        %Set average reference EEG data equal to record
        eval(['record = Session_',num2str(raw_files),'_EEG_avgref;'])

        %Set bad samples to rej_samples
        eval(['rej_samples = Session_',num2str(raw_files),'_bad_samples;'])

        %Load event time data and trial identifiers
        cd(events_dir)
        load(['Raw_file_', num2str(raw_files),'_DIN_marker_times.mat']);
        load(['Raw_file_', num2str(raw_files),'_trial_identifiers.mat']);

        %Remove the last trial from 600 because missing DINs 
        if isequal(ID, '600') && isequal(raw_files, 2) && isequal(relevant_location, 'Center Relevant')
            
            first_stim_identifier(72) = [];
            
        end
        
        %% Select the task relevant and irrelevant stimuli

        %NOTE: DIN marker times are defined as either faces, questions,
        %button presses, and trial start times. This next step of the
        %analysis is to further subdivide the face events as either
        %relevant faces (reported) or irrelevant faces (not reported).
        %Then those events can be defined among opacity and perception
        %metrics.

        %Initialize variables
        DIN_relevant_faces = [];
        DIN_irrelevant_faces = [];

        %Loop over rows
        for row = 1:size(first_stim_identifier,1)

            %Quadrant relevant
            if isequal(relevant_location, 'Quadrant Relevant')

                %Center First Trial
                if strcmp(first_stim_identifier{row}, 'Face Center')

                    %Relevant face - Get the second face time in the trial 
                    DIN_relevant_faces = [DIN_relevant_faces; DIN_faces(row*2)];

                    %Irrelevant face - Get the first face time in the trial
                    DIN_irrelevant_faces = [DIN_irrelevant_faces; DIN_faces(row*2-1)];

                %Quadrant First Trial
                elseif strcmp(first_stim_identifier{row}, 'Face Quadrant')

                    %Relevant face - Get the first face time in the trial
                    DIN_relevant_faces = [DIN_relevant_faces; DIN_faces(row*2-1)];

                    %Irrelevant face - Get the first face time in the trial
                    DIN_irrelevant_faces = [DIN_irrelevant_faces; DIN_faces(row*2)];

                end

            %Center relevant    
            elseif isequal(relevant_location, 'Center Relevant')

                %Center First
                if strcmp(first_stim_identifier{row}, 'Face Quadrant')

                    %Relevant face - Get the second face time in the trial
                    DIN_relevant_faces = [DIN_relevant_faces; DIN_faces(row*2)];

                    %Irrelevant face - Get the first face time in the trial
                    DIN_irrelevant_faces = [DIN_irrelevant_faces; DIN_faces(row*2-1)];

                %Quadrant First
                elseif strcmp(first_stim_identifier{row}, 'Face Center')

                    %Relevant face - Get the first face time in the trial
                    DIN_relevant_faces = [DIN_relevant_faces; DIN_faces(row*2-1)];

                    %Irrelevant face - Get the first face time in the trial
                    DIN_irrelevant_faces = [DIN_irrelevant_faces; DIN_faces(row*2)];

                end

            end

        end

        %Loop through the epoch event types
        for event = 1:size(epoch_types,2)

            %Create index variable for each event type
            if strcmp(epoch_types{event},'relevant_faces') == 1
                IDX = DIN_relevant_faces;

            elseif strcmp(epoch_types{event},'irrelevant_faces') == 1
                IDX = DIN_irrelevant_faces;

            elseif strcmp(epoch_types{event},'questions') == 1
                IDX = DIN_questions;

            elseif strcmp(epoch_types{event},'button_presses') == 1
                IDX = DIN_buttons;
                
            elseif strcmp(epoch_types{event},'relevant_faces_post_jitter') == 1
                IDX = DIN_relevant_faces + 3000; %Add 3s from stimulus onset
            
            elseif strcmp(epoch_types{event},'irrelevant_faces_post_jitter') == 1
                IDX = DIN_irrelevant_faces + 3000; %Add 3s from stimulus onset
            
            end

            %% Cut Epochs of Specified Size and Reject Epochs with too many bad samples

            disp(['Cutting ', epoch_types{event},' ', num2str(half_epoch_duration*2), ' epochs'])

            %Create an empty epochs matrix
            epochs = [];
            clean_samples_epochs = [];
            bad_epochs_idx = [];

            % Loop through the number of IDX values (i.e., the DIN markers for a particular event type)
            for i = 1:size(IDX,1)

                % Set epoch to the EEG data that is half_epoch_size before and after
                % IDX(i) or a particular face event
                trial = record(:,IDX(i)-sampling_rate*half_epoch_duration:IDX(i)+sampling_rate*half_epoch_duration); 

                %Cut reject sample epochs (0 = bad sample; 1 = good sample)
                clean_samples_trial = rej_samples(1,IDX(i)-sampling_rate*half_epoch_duration:IDX(i)+sampling_rate*half_epoch_duration);

                %Aggregate rej_epoch into time x trial matrix
                clean_samples_epochs(:,i) = clean_samples_trial;

                %Find bad samples and reject epochs with too many (>25%
                %samples) during the bad window (-200 to +500ms)
                if length(find(clean_samples_trial(bad_window)<1))/length(bad_window) < bad_sample_prop

                    %Enter 0 for a good epoch
                    bad_epochs_idx(i,1) = 0;

                elseif length(find(clean_samples_trial(bad_window)<1))/length(bad_window) >= bad_sample_prop

                    %Enter 1 for a bad epoch
                    bad_epochs_idx(i,1) = 1;

                end

                %Normalized voltage to the average of each channel signal within the baseline period of the epoch.  
                %(1) Find the nanmean of epoch across time for 1 to length of fs (or the baseline period of the epoch prior to the face)
                %(2) Replicate the mean value for each channel to the length of epoch
                %(3) Subtract each value in epoch with the corresponding
                %channel mean (trial level baselining - baseline is the entire first 2000ms pre stimulus)
                epochs(:,:,i) = (trial-repmat(nanmean(trial(:,1:sampling_rate*half_epoch_duration),2),[1,size(trial,2)]));

            end

            %Rename epochs, rej_epochs_idx, and rej_epochs variable to
            %their Raw/session file
            eval(['Raw_file_',num2str(raw_files),'_',epoch_types{event}, '_epochs = epochs;'])
            eval(['Raw_file_',num2str(raw_files),'_',epoch_types{event}, '_rej_epochs = clean_samples_epochs;'])
            eval(['Raw_file_',num2str(raw_files),'_',epoch_types{event}, '_bad_epochs_idx = bad_epochs_idx;'])

         end

    end

    %% Combine Epochs Across Raw Files

    disp('Combining epochs across raw session files')

    %Loop over epochs types
    for type = 1:size(epoch_types,2)

        %Initialize variable
        all_epochs = [];
        all_clean_samples_epochs = [];
        all_bad_epochs_idx = [];

        %Select event type
        event_type = epoch_types{type};

        %Loop over sessions and combine epochs
        for raw_files = 1:size(EEG_data,1)

            %Rename epochs to generic variable name
            session_epochs = eval(['Raw_file_',num2str(raw_files),'_',event_type, '_epochs']); 
            session_clean_samples_epochs = eval(['Raw_file_',num2str(raw_files),'_',event_type, '_rej_epochs']);
            session_bad_epochs_idx = eval(['Raw_file_',num2str(raw_files),'_',event_type, '_bad_epochs_idx']);

            %Combine epochs across sessions
            all_epochs = cat(3,all_epochs, session_epochs); %channels x time x trials
            all_clean_samples_epochs = cat(2,all_clean_samples_epochs, session_clean_samples_epochs); %time x trials
            all_bad_epochs_idx = cat(1,all_bad_epochs_idx, session_bad_epochs_idx); %trials

        end

        %Rename combined variables
        eval(['All_',event_type,'_epochs = all_epochs;'])
        eval(['All_',event_type,'_epochs_clean_samples_mask = all_clean_samples_epochs;'])
        eval(['All_',event_type,'_bad_epochs_idx = all_bad_epochs_idx;'])

    end

    %Save combined main event epochs
    cd(preprocessed_dir)
    save(['Avg_ref_main_events_',num2str(half_epoch_duration*2),'s_epochs.mat'], 'All*epochs', '-append');
    
    %Save combined index of bad epoch index and clean sample mask epochs
    cd(bad_channel_dir)
    save('Main_events_bad_epochs_samples.mat', 'All*epochs_clean_samples_mask', 'All*bad_epochs_idx', '-append');

end