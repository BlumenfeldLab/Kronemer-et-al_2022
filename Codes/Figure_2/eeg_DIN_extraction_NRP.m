%% EEG DIN Marker Extraction

%Finds the DIN marker times from the EEG data files and correct the DIN
%time by 42ms based on photodiode testing.

%Written by: Sharif I. Kronemer
%Date: 2/7/2021

function eeg_DIN_extraction_NRP(ID, Beh_dir, events_dir, EEG_dir, relevant_location)

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

        %% Find DIN Marker Times Generated in Net Station
        
        %Updating event times in EEG structure that are lost in certain
        %participants
        if isequal(ID,'567') && isequal(raw_files,2) && isequal(relevant_location, 'Center Relevant')
            
            EEG.event = [EEG.event(170:177),EEG.event];
            
        elseif isequal(ID,'599') && isequal(raw_files,2) && isequal(relevant_location, 'Quadrant Relevant')
            
            EEG.event = [EEG.event(170:177),EEG.event];
                        
        end  

        %% DIN Faces Events (Fac2)

        DIN_faces = [];  

        %Loop through the number of DIN markers
        for row = 7:size(EEG.event,2)-1 %Note: skipping over initial DIN markers 

            %Find first face stimulus in a trial
            if strcmp(EEG.event(row).type,'Fac2') == 1 && strcmp(EEG.event(row-1).type,'Trl1') == 1 && strcmp(EEG.event(row-2).type,'Pre8') == 1 %Identifies the DIN marker sequence that corresponds with face stimuli 

               %Saves the DIN time as "event"
               event = [EEG.event(row).latency]; 
               DIN_faces = [DIN_faces; event]; 

            %Find second face stimulus in a trial   
            elseif strcmp(EEG.event(row).type,'Fac2') == 1 && strcmp(EEG.event(row+1).type,'Que4') == 1 && strcmp(EEG.event(row+2).type,'Pre8') == 1 %Identifies the DIN marker sequence that corresponds with face stimuli 

               %Saves the DIN time as "event"
               event = [EEG.event(row).latency]; 
               DIN_faces = [DIN_faces; event]; 

            %Find the first stimulus in the first trial of a run
            elseif strcmp(EEG.event(row).type,'Fac2') == 1 && strcmp(EEG.event(row-1).type,'Trl1') == 1 && strcmp(EEG.event(row-6).type,'Trl1') == 1 %Identifies the DIN marker sequence that corresponds with face stimuli 

               %Saves the DIN time as "event"
               event = [EEG.event(row).latency]; 
               DIN_faces = [DIN_faces; event];    

            else     
                continue

            end

        end
        
        %If subjects have usual trial or session counts skip from
        %determining if they have the correct number of trials
        if isequal(ID, '608') && isequal(relevant_location, 'Quadrant Relevant') ||...
            isequal(ID, '607') && isequal(relevant_location, 'Quadrant Relevant') ||...
              isequal(ID, '652') && isequal(relevant_location, 'Center Relevant') ||...
                isequal(ID, '623') && isequal(relevant_location, 'Quadrant Relevant') ||...
                    isequal(ID, '600') && isequal(relevant_location, 'Center Relevant')
                
            disp('Skipping face trial count confirmation')
            
        else
            
            %Confirm the number of faces is correct
            if isequal(raw_files, 1) && ~isequal(size(DIN_faces,1), 96)

                error('Number of faces incorrect in Raw file 1')

            elseif isequal(raw_files, 2) && ~isequal(size(DIN_faces,1), 144)

                error('Number of faces incorrect in Raw file 2')

            end
            
        end

        %% DIN Questions Events (Que4)

        DIN_questions = [];

        %Loop through the number of DIN markers
        for row = 1:size(EEG.event,2)  

            %Confirms that the correct question events are selected and not thoses corresponding to trial onset and offset where all pulse types are sent
            if strcmp(EEG.event(row).type,'Que4') == 1 && strcmp(EEG.event(row-1).type,'Fac2') == 1 && strcmp(EEG.event(row+1).type,'Pre8') == 1 %Identifies the DIN marker sequence that corresponds with the perception question

                %Saves the DIN time as "event"
                event = [EEG.event(row).latency]; 
                DIN_questions = [DIN_questions; event]; 

            elseif strcmp(EEG.event(row).type,'Que4') == 1 && strcmp(EEG.event(row-1).type,'Pre8') == 1 && strcmp(EEG.event(row+1).type,'Pre8') == 1 %Identifies the DIN marker sequence that corresponds with the location question

                %Saves the DIN time as "event"
                event = [EEG.event(row).latency]; 
                DIN_questions = [DIN_questions; event]; 

            else 
                continue

            end

        end

        %If subjects have usual trial or session counts skip from
        %determining if they have the correct number of trials
        if isequal(ID, '608') && isequal(relevant_location, 'Quadrant Relevant') ||...
            isequal(ID, '607') && isequal(relevant_location, 'Quadrant Relevant') ||...
             isequal(ID, '652') && isequal(relevant_location, 'Center Relevant') ||...
               isequal(ID, '623') && isequal(relevant_location, 'Quadrant Relevant') ||...
                isequal(ID, '600') && isequal(relevant_location, 'Center Relevant')
    
            disp('Skipping question trial count confirmation')
            
        else
            
            %Confirm the number of faces is correct
            if isequal(raw_files, 1) && ~isequal(size(DIN_questions,1), 96)

                error('Number of questions incorrect in Raw file 1')

            elseif isequal(raw_files, 2) && ~isequal(size(DIN_questions,1), 144)

                error('Number of questions incorrect in Raw file 2')

            end

        end
        
        %% DIN Button Press Events (Pre8)

        DIN_buttons = [];

        %Loop through the number of DIN markers
        for row = 1:size(EEG.event,2)

            %Skip first row for 458DW
            if isequal(row,1) && isequal(raw_files, 2) && isequal('458DW',ID)

                continue

            end

            %If the DIN marker is not the last DIN marker
            if row < size(EEG.event,2)

                %Confirms that the correct button press events are selected and not thoses corresponding to trial onset and offset where all pulse types are sent
                if strcmp(EEG.event(row).type,'Pre8') == 1 && strcmp(EEG.event(row-1).type,'Que4') == 1 && strcmp(EEG.event(row+1).type,'Que4')== 1

                   %Saves the DIN time as "event"     
                   event = [EEG.event(row).latency];
                   DIN_buttons = [DIN_buttons; event];

                elseif strcmp(EEG.event(row).type,'Pre8') == 1 && strcmp(EEG.event(row-1).type,'Que4')== 1 && strcmp(EEG.event(row+1).type,'Trl1')== 1

                   %Saves the DIN time as "event" 
                   event = [EEG.event(row).latency];
                   DIN_buttons = [DIN_buttons; event];

                else 
                    continue

                end  

            %If the DIN marker is the last DIN marker   
            elseif row == size(EEG.event,2)  

               if strcmp(EEG.event(row).type,'Pre8') == 1 && strcmp(EEG.event(row-1).type,'Que4')== 1 && strcmp(EEG.event(row-2).type,'Pre8')== 1

                   %Saves the DIN time as "event"
                   event = [EEG.event(row).latency];
                   DIN_buttons = [DIN_buttons; event];

               end

            end

        end

        %If subjects have usual trial or session counts skip from
        %determining if they have the correct number of trials
        if isequal(ID, '608') && isequal(relevant_location, 'Quadrant Relevant') ||...
            isequal(ID, '607') && isequal(relevant_location, 'Quadrant Relevant') ||...
             isequal(ID, '652') && isequal(relevant_location, 'Center Relevant') ||...
               isequal(ID, '623') && isequal(relevant_location, 'Quadrant Relevant') ||...
                isequal(ID, '600') && isequal(relevant_location, 'Center Relevant')
        
            disp('Skipping button press trial count confirmation')
            
        else
            
            %Confirm the number of faces is correct
            if isequal(raw_files, 1) && ~isequal(size(DIN_buttons,1), 96)

                error('Number of buttons incorrect in Raw file 1')

            elseif isequal(raw_files, 2) && ~isequal(size(DIN_buttons,1), 144)

                error('Number of buttons incorrect in Raw file 2')

            end

        end
        
        %% DIN Trials (Trl1)

        DIN_trials = [];

        %Loop through rows 
        for row = 6:size(EEG.event,2)-1 %Note: starting from 6th row to skip over initial DIN markers  

            %Confirms that the correct question events are selected and not thoses corresponding to trial onset and offset where all pulse types are sent
            if strcmp(EEG.event(row).type,'Trl1') == 1 && strcmp(EEG.event(row-1).type,'Pre8') == 1 && strcmp(EEG.event(row+1).type,'Fac2') == 1 %Identifies the DIN marker sequence that corresponds with the perception question

                %Saves the DIN time as "event"
                event = [EEG.event(row).latency]; 
                DIN_trials = [DIN_trials; event]; 

            %Find the first trial in a run    
            elseif strcmp(EEG.event(row).type,'Trl1') == 1 && strcmp(EEG.event(row-5).type,'Trl1') == 1 && strcmp(EEG.event(row+1).type,'Fac2') == 1 %Identifies the DIN marker sequence that corresponds with the perception question

                %Saves the DIN time as "event"
                event = [EEG.event(row).latency]; 
                DIN_trials = [DIN_trials; event]; 

            else 
                continue

            end

        end
        
        %If subjects have usual trial or session counts skip from
        %determining if they have the correct number of trials
        if isequal(ID, '608') && isequal(relevant_location, 'Quadrant Relevant') ||...
            isequal(ID, '607') && isequal(relevant_location, 'Quadrant Relevant') ||...
             isequal(ID, '652') && isequal(relevant_location, 'Center Relevant') ||...
               isequal(ID, '623') && isequal(relevant_location, 'Quadrant Relevant')
        
            disp('Skipping trial count confirmation')
            
        else
            
            %Confirm the number of faces is correct
            if isequal(raw_files, 1) && ~isequal(size(DIN_trials,1), 48)

                error('Number of trials incorrect in Raw file 1')

            elseif isequal(raw_files, 2) && ~isequal(size(DIN_trials,1), 72)

                error('Number of trials incorrect in Raw file 2')

            end

        end
        
        %% *** Correct TTL and Photodiode Times ****

        %NOTE: There is a lag between the TTL pulse up-phase time and the
        %photodiode up-phase which can range between 40-60ms. Adding a
        %constant value to the face events will correct for this lag. 

        % Add 42 ms to face time
        DIN_faces = DIN_faces + 42; 

        %% Save the DIN Marker Times

        cd(events_dir)
        save(['Raw_file_', num2str(raw_files),'_DIN_marker_times.mat'], 'DIN_faces', 'DIN_questions', 'DIN_buttons', 'DIN_trials')

    end

end