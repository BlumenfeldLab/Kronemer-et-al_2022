%% Cut Specific Epochs from the Main Event Epochs

%This function takes the ICA component removed epochs and cuts specfic
%epoch types using the behavioral identifier information, including CP,
%CnP, and irrelevant stimuli by opacity type. This epochs are saved in the
%epochs directory. Bad epochs identified earlier for each main event type
%are excluded. 

%Written by: Sharif Kronemer
%Date: 2/6/2021
%Modified: 4/29/2021

function cut_eeg_specific_epochs_report_paradigm(EEG_data, events_dir, epochs_dir, bad_channel_dir, half_epoch_duration)
    
    %Loop through the different epoch lengths
    epoch_length = [num2str(half_epoch_duration*2),'s'];
    
    %Load general epochs
    cd(epochs_dir)
    load('face_epochs_preprocessed_components_removed.mat');
    
    %Load bad epochs index
    cd(bad_channel_dir)
    load('Main_events_bad_epochs_samples.mat')
    
    %Initialize variables
    all_perception_identifier = [];
    all_delay_identifier = [];

    %Loop over raw files and combine the identifiers
    for raw_files = 1:size(EEG_data,2)

        %Load event time data
        cd(events_dir)
        load(['Raw_file_', num2str(raw_files+1),'_trial_identifiers.mat']);
        
        %Combine identifiers 
        all_perception_identifier = [all_perception_identifier; perception_identifier];
        all_delay_identifier = [all_delay_identifier; delay_identifier];
        
    end
    
    %% Categorize all the face event types by perception 

    disp('Cutting specific epochs')

    %Create empty matrices for each perception type and post-stim times
    
    %Relevant face variables
    CP_15s_all_epochs_idx = []; %Confirmed perceived threshold
    CP_15s_bad_epochs_idx = []; %Confirmed perceived threshold
    
    CP_1s_all_epochs_idx = []; %Confirmed perceived opaque
    CP_1s_bad_epochs_idx = []; %Confirmed perceived opaque

    CnP_15s_all_epochs_idx = []; %Confirmed not perceived threshold 
    CnP_15s_bad_epochs_idx = []; %Confirmed not perceived threshold 

    CnP_1s_all_epochs_idx = []; %Confirmed not perceived opaque
    CnP_1s_bad_epochs_idx = []; %Confirmed not perceived opaque

    CG_all_epochs_idx = []; %Correct guess (participant said she did not see stimulus but accurately locates stimulus)
    CG_bad_epochs_idx = []; %Correct guess (participant said she did not see stimulus but accurately locates stimulus)   
    
    FPer_all_epochs_idx = []; %False perception (participant see stimulus but locates inaccurately) 
    FPer_bad_epochs_idx = []; %False perception (participant see stimulus but locates inaccurately) 

    %Loop through the perception_identifier variable and categorize by type
    %Creates matrix of the row number with certain event types. 
    for row = 1:size(all_perception_identifier, 1)
        
        % Confirmed perceived threshold opacity trials
        if strcmp('Confirmed perceived', all_perception_identifier(row)) == 1 && all_delay_identifier(row) == 15 
            
            CP_15s_all_epochs_idx = [CP_15s_all_epochs_idx; row]; %Sets variable to the row number corresponding to the face events in sequences (i.e., row 3 of the identifier variables = face event 3)

            %Determine if bad relevant trial
            if isequal(All_faces_bad_epochs_idx(row),1)
                
                %Store index of bad irrelevant relative to this epoch type
                CP_15s_bad_epochs_idx = [CP_15s_bad_epochs_idx; size(CP_15s_all_epochs_idx,1)];
        
            end            
            
        % Confirmed perceived opaque opacity trails    
        elseif strcmp('Confirmed perceived', all_perception_identifier(row)) == 1 && all_delay_identifier(row) == 1  
            
            CP_1s_all_epochs_idx = [CP_1s_all_epochs_idx; row]; %Sets variable to the row number corresponding to the face events in sequences (i.e., row 3 of the identifier variables = face event 3)

            %Determine if bad relevant trial
            if isequal(All_faces_bad_epochs_idx(row),1)
                
                %Store index of bad irrelevant relative to this epoch type
                CP_1s_bad_epochs_idx = [CP_1s_bad_epochs_idx; size(CP_1s_all_epochs_idx,1)];
        
            end     
            
        % Confirmed not perceived trials
        elseif strcmp('Confirmed not perceived', all_perception_identifier(row)) == 1 && all_delay_identifier(row) == 15 
            
            CnP_15s_all_epochs_idx = [CnP_15s_all_epochs_idx; row]; %Sets variable to the row number corresponding to the face events in sequences (i.e., row 3 of the identifier variables = face event 3)

            %Determine if bad relevant trial
            if isequal(All_faces_bad_epochs_idx(row),1)
                
                %Store index of bad irrelevant relative to this epoch type
                CnP_15s_bad_epochs_idx = [CnP_15s_bad_epochs_idx; size(CnP_15s_all_epochs_idx,1)];
        
            end     
            
        % Confirmed not perceived trials
        elseif strcmp('Confirmed not perceived', all_perception_identifier(row)) == 1 && all_delay_identifier(row) == 1  
            
            CnP_1s_all_epochs_idx = [CnP_1s_all_epochs_idx; row]; %Sets variable to the row number corresponding to the face events in sequences (i.e., row 3 of the identifier variables = face event 3)    

            %Determine if bad relevant trial
            if isequal(All_faces_bad_epochs_idx(row),1)
                
                %Store index of bad irrelevant relative to this epoch type
                CnP_1s_bad_epochs_idx = [CnP_1s_bad_epochs_idx; size(CnP_1s_all_epochs_idx,1)];
        
            end     
            
        % Correct guess trials
        elseif strcmp('Correct guess', all_perception_identifier(row)) == 1
            
            CG_all_epochs_idx = [CG_all_epochs_idx; row]; %Sets variable to the row number corresponding to the face events in sequences (i.e., row 3 of the identifier variables = face event 3)

            %Determine if bad relevant trial
            if isequal(All_faces_bad_epochs_idx(row),1)
                
                %Store index of bad irrelevant relative to this epoch type
                CG_bad_epochs_idx = [CG_bad_epochs_idx; size(CG_all_epochs_idx,1)];
        
            end     
            
        % False perception trials    
        elseif strcmp('False perception', all_perception_identifier(row)) == 1 
            
            FPer_all_epochs_idx = [FPer_all_epochs_idx; row]; %Sets variable to the row number corresponding to the face events in sequences (i.e., row 3 of the identifier variables = face event 3)
            
            %Determine if bad relevant trial
            if isequal(All_faces_bad_epochs_idx(row),1)
                
                %Store index of bad irrelevant relative to this epoch type
                FPer_bad_epochs_idx = [FPer_bad_epochs_idx; size(FPer_all_epochs_idx,1)];
        
            end    
            
        end

    end

    %% Cut epochs for each face event identifier (e.g., CP, CnP, CG, etc.)
    
    %% Confirmed Perceived 1s

    %Create empty vector 
    CP_1s_epochs = [];

    %Populate vector with data
    CP_1s_epochs = face_epochs_components_removed_data(:,:,CP_1s_all_epochs_idx);

    %Save
    cd(epochs_dir)
    save(['CP_1s_',epoch_length,'_epochs.mat'], 'CP_1s_epochs', 'CP_1s_all_epochs_idx') 

    %% Confirmed Perceived 15s

    %Create empty vector 
    CP_15s_epochs = [];

    %Populate vector with data
    CP_15s_epochs = face_epochs_components_removed_data(:,:,CP_15s_all_epochs_idx);

    %Save
    cd(epochs_dir)
    save(['CP_15s_',epoch_length,'_epochs.mat'], 'CP_15s_epochs', 'CP_15s_all_epochs_idx') 

    %% Confirmed not Perceived 1s

    %Create empty vector 
    CnP_1s_epochs = [];

    %Populate vector with data
    CnP_1s_epochs = face_epochs_components_removed_data(:,:,CnP_1s_all_epochs_idx);

    %Save
    cd(epochs_dir)
    save(['CnP_1s_',epoch_length,'_epochs.mat'], 'CnP_1s_epochs', 'CnP_1s_all_epochs_idx') 

    %% Confirmed not Perceived 15s

    %Create empty vector 
    CnP_15s_epochs = [];

    %Populate vector with data
    CnP_15s_epochs = face_epochs_components_removed_data(:,:,CnP_15s_all_epochs_idx);

    %Save
    cd(epochs_dir)
    save(['CnP_15s_',epoch_length,'_epochs.mat'], 'CnP_15s_epochs', 'CnP_15s_all_epochs_idx') 

   %% Correct Guess

    %Create empty vector 
    CG_epochs = [];

    %Populate vector with data
    CG_epochs = face_epochs_components_removed_data(:,:,CG_all_epochs_idx);

    %Save
    cd(epochs_dir)
    save(['CG_',epoch_length,'_epochs.mat'], 'CG_epochs', 'CG_all_epochs_idx') 

    %% False Perception

    %Create empty vector 
    FPer_epochs = [];

    %Populate vector with data
    FPer_epochs = face_epochs_components_removed_data(:,:,FPer_all_epochs_idx);

    %Save
    cd(epochs_dir)
    save(['FPer_',epoch_length,'_epochs.mat'], 'FPer_epochs', 'FPer_all_epochs_idx') 
    
%{
    %% Questions
    
    %Rename data and remove bad question epochs
    questions_epochs = questions_epochs_components_removed;
    %questions_epochs = questions_epochs_components_removed(:,:,~All_questions_bad_epochs_idx);
    
    %Rename bad epochs index for consistency with other variable types
    questions_bad_epochs_idx = find(All_questions_bad_epochs_idx);
    
    %Save
    cd(epochs_dir)
    save(['questions_',epoch_length,'_epochs.mat'], 'questions_epochs') 
    
    %% Button Presses
    
    %Rename data and remove bad button press epochs
    button_presses_epochs = button_presses_epochs_components_removed;
    %button_presses_epochs = button_presses_epochs_components_removed(:,:,~All_button_presses_bad_epochs_idx);
    
    %Rename bad epochs index for consistency with other variable types
    button_presses_bad_epochs_idx = find(All_button_presses_bad_epochs_idx);
    
    %Save
    cd(epochs_dir)
    save(['button_presses_',epoch_length,'_epochs.mat'], 'button_presses_epochs') 
%}    
    
    %% Save bad epochs index
    
    cd(events_dir)
    save('Bad_trials_index.mat', 'CP_1s_bad_epochs_idx', 'CP_15s_bad_epochs_idx', 'CnP_1s_bad_epochs_idx', ...
        'CnP_15s_bad_epochs_idx')
    
end