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
%Modified: 4/29/2019

%% Loop through all the raw files available for testing session
 
function eeg_DIN_extraction_report_paradigm(ID, EEG_dir, EEG_data, EEG_folders, events_dir, run_location)
   
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
            
            EEG = pop_readegi(fullfile(EEG_dir, EEG_folders(raw_files).name, EEG_data{raw_files}));
            
        else
            
            EEG = pop_readegi(fullfile(EEG_dir, EEG_folders(raw_files+1).name, EEG_data{raw_files}));
            
        end
        
        %Error of data struct is empty
        if isempty(EEG) 
            
            error('EEG structure is empty') 
       
        end 

        %% Find DIN (EEG ditial inpute markers) Marker Times Generated in Net Station

        %% DIN Faces Events (Fac2)
        
        DIN_faces = [];

            %Loop through the number of DIN markers
            for row = 1:size(EEG.event,2) 

                %Confirms that the correct face events are selected and not thoses corresponding to trial onset and offset where all pulse types are sent
                if strcmp(EEG.event(row).type,'Fac2') == 1 && strcmp(EEG.event(row-1).type,'Trl1') == 1 && strcmp(EEG.event(row+1).type,'Que4') == 1 %Identifies the DIN marker sequence that corresponds with face stimuli 

                   %Saves the DIN time as "event"
                   event = [EEG.event(row).latency]; 

                   %Adds the "event" value to EEG_faces
                   DIN_faces(end+1,1) = event; 

                else     
                    continue
                end

            end 

        %% DIN Questions Events (Que4)
        
        DIN_questions = [];

            for row = 1:size(EEG.event,2)  

                %Confirms that the correct question events are selected and not thoses corresponding to trial onset and offset where all pulse types are sent
                if strcmp(EEG.event(row).type,'Que4') == 1 && strcmp(EEG.event(row-1).type,'Fac2') == 1 && strcmp(EEG.event(row+1).type,'Pre8') == 1 %Identifies the DIN marker sequence that corresponds with the perception question

                    %Saves the DIN time as "event"
                    event = [EEG.event(row).latency]; 

                    %Adds the "event" value to EEG_questions
                    DIN_questions(end+1,1) = event; 

                elseif strcmp(EEG.event(row).type,'Que4') == 1 && strcmp(EEG.event(row-1).type,'Pre8') == 1 && strcmp(EEG.event(row+1).type,'Pre8') == 1 %Identifies the DIN marker sequence that corresponds with the location question

                    %Saves the DIN time as "event"
                    event = [EEG.event(row).latency]; 

                    %Adds the "event" value to EEG_questions
                    DIN_questions(end+1,1) = event; 

                else 
                    continue
                end

            end

        %% DIN Button Press Events (Pre8)
        
        DIN_buttons = [];

            for row = 1:size(EEG.event,2)

                %If the DIN marker is not the last DIN marker
                if row < size(EEG.event,2)

                    %Confirms that the correct button press events are selected and not thoses corresponding to trial onset and offset where all pulse types are sent
                    if strcmp(EEG.event(row).type,'Pre8') == 1 && strcmp(EEG.event(row-1).type,'Que4') == 1 && strcmp(EEG.event(row+1).type,'Que4')== 1

                       event = [EEG.event(row).latency];
                       DIN_buttons(end+1,1) = event;

                    elseif strcmp(EEG.event(row).type,'Pre8') == 1 && strcmp(EEG.event(row-1).type,'Que4')== 1 && strcmp(EEG.event(row+1).type,'Trl1')== 1

                       event = [EEG.event(row).latency];
                       DIN_buttons(end+1,1) = event;

                    else 
                        continue
                    end  

                %If the DIN marker is the last DIN marker   
                elseif row == size(EEG.event,2)  

                   if strcmp(EEG.event(row).type,'Pre8') == 1 && strcmp(EEG.event(row-1).type,'Que4')== 1 && strcmp(EEG.event(row-2).type,'Pre8')== 1

                       event = [EEG.event(row).latency];
                       DIN_buttons(end+1,1) = event;

                   end

                end

            end

        %% Subject specific correction - Remove the last face and last two questions and button presses DIN markers for 260DY due to a recording error    

        if strcmp(ID, '260DY')

           if size(DIN_faces,1) == 33   
              DIN_faces(33) =[];   
           end

           if size(DIN_buttons,1) == 66
               DIN_buttons(65:66) = [];
           end

           if size(DIN_questions,1) == 66
               DIN_questions(65:66) = []; 
           end

        end

        %% Subject specific correction - Remove the last face event and corresponding button events from trial 32 of run 1 for 243EG

        if strcmp(ID, '243EG')

            if isequal(raw_files, 1) == 1

               if size(DIN_faces,1) == 64   
                  DIN_faces(32) =[];   
               end

               if size(DIN_buttons,1) == 128
                   DIN_buttons(63:64) = [];
               end

               if size(DIN_questions,1) == 128
                   DIN_questions(63:64) = []; 
               end

            end

        end

        %% ****Correct for Arduino TTL Pulse Time to EEG System****

        %NOTE: TTL pulses sent to the EEG system are delayed by the duration of
        %the TTL pulse because pulses are first sent to the EyeLink
        %(pupillometry) system. Therefore, the EEG amplifier is receiving TTL
        %pulses that are delayed by up to 300ms. This portion of the code
        %corrects for this lag by subtracting this delay from the DIN marker
        %times. 

        % Face events TTL pulses must be corrected by 151ms (face TTL pulse =
        % 150ms)
        DIN_faces = DIN_faces - 151;

        % Question events TTL pulses must be corrected by 101ms (question TTL
        % pulse = 100ms)
        DIN_questions = DIN_questions - 101;

        % Button press events TTL pulses must be corrected by 26ms (question
        % TTL pulse = 26ms)
        DIN_buttons = DIN_buttons - 26;

        %% *** Correct TTL and Photodiode Times ****

        %NOTE: There is a lag between the TTL pulse up-phase time and the
        %photodiode up-phase which can range between 40-60ms. Adding a
        %constant value to the face events will correct for this lag. 

        % Add 42 ms to face time
        DIN_faces = DIN_faces + 42; 

        %% Save the DIN Marker Times

        cd(events_dir)
        save(['Raw_file_', num2str(raw_files+1),'_DIN_marker_times.mat'], 'DIN_faces', 'DIN_questions', 'DIN_buttons')
        
        %% Compare PsychoPy and DIN Marker Times and Confirm Events Correspond Using the Function verify_psychopy_ttl_sEEG_11_9_2017
 
        %{    
        % For each event time two figures will be created and saved showing the
        % relationship between the DIN times and Psychopy. The experimenter
        % should review each of these figures to confirm that the timing
        % characteristics of the are ideal. For example, the clock drift plots
        % the the difference in time between the PsychoPy and Net Station
        % computers should increase linearly (non-linear changes are
        % unacceptable). 

        %% Face onset pulses

        %Compares the face DIN markers versus Psychopy indices 
        event  = 'Faces';
        [face_mismatches, face_explanation] = verify_psychopy_ttl_sEEG_11_9_2017(EEG_faces, psychopy_faces, events_dir, event, raw_files,ID);

        %% Question onset pulses

        %Compares the question DIN markers versus Psychopy indices 
        event = 'Questions';
        [question_mismatches, question_explanation] = verify_psychopy_ttl_sEEG_11_9_2017(EEG_questions, reshape(psychopy_questions,[1 numel(psychopy_questions)]),events_dir,event,raw_files,ID);

        %% Button press pulses

        %Compares the button DIN markers versus Psychopy indices 
        event = 'Button presses';
        [button_mismatches, button_explanation] = verify_psychopy_ttl_sEEG_11_9_2017(EEG_buttons, reshape(psychopy_buttons,[1 numel(psychopy_buttons)]),events_dir,event,raw_files,ID);
        %}
        
    end