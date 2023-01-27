%% Identify CP and CnP Trials from Behavioral Log File

%The goal of this code is to load the behavioral log file and identify CP
%and CnP trials from the events registered with the log file. This code was
%particularly designed for the visual no report paradigm version 1.12. 

%Note: This code has been modified to allow managing multiple log files
%that are created when the behavioral task is restarted during a testing
%session.

%Written by: Sharif I. Kronemer
%Date: 7/15/2019
%Modified: 1/5/2021

function Find_CP_CnP_Events_From_Log_File_Function(Yes, beh_folder, savefolder)

    %Load the behavioral log file
    cd(beh_folder)
    
    %Find the log files in behavioral directory
    log_files = dir('*.log'); 
    
    %Count the number of log files
    num_log_files = size(log_files,1); 

    %% Initialize variables

    %Stimuli time
    all_relevant_face_times = [];
    all_irrelevant_face_times = [];

    %Stimuli position
    relevant_face_position = [];
    irrelevant_face_position = [];

    %Stimuli opacity 
    relevant_face_opacity = [];
    irrelevant_face_opacity = [];
    current_relevant_face_pos = ['empty'];

    %Keypress
    all_per_ques_keypresses = [];
    all_loc_ques_keypresses = [];

    %CP threshold run variables
    CP_thres_face_times_run_1 = [];
    CP_thres_face_times_run_2 = [];
    CP_thres_face_times_run_3 = [];
    CP_thres_face_times_run_4 = [];
    CP_thres_face_times_run_5 = [];

    %CnP threshold run variables
    CnP_thres_face_times_run_1 = [];
    CnP_thres_face_times_run_2 = [];
    CnP_thres_face_times_run_3 = [];
    CnP_thres_face_times_run_4 = [];
    CnP_thres_face_times_run_5 = [];

    %CP opaque run variables
    CP_opaque_face_times_run_1 = [];
    CP_opaque_face_times_run_2 = [];
    CP_opaque_face_times_run_3 = [];
    CP_opaque_face_times_run_4 = [];
    CP_opaque_face_times_run_5 = [];

    %CnP opaque run variables
    CnP_opaque_face_times_run_1 = [];
    CnP_opaque_face_times_run_2 = [];
    CnP_opaque_face_times_run_3 = [];
    CnP_opaque_face_times_run_4 = [];
    CnP_opaque_face_times_run_5 = [];

    %% Loop over log files
    for file = 1:num_log_files
        
        %Load current log data file
        DTA = importdata(log_files(file).name);

        %% Find the trial onset times

        %Initialize variable
        trial_start_times = [];

        %Loop over rows of log file
        for row = 1:length(DTA)

            %Search for string indicating start of a trial
            if  any(~cellfun('isempty',strfind(DTA(row),'Starting trial'))) 

               %Create index of trial start rows
               trial_start_times = cat(1, trial_start_times, row); 

            end

        end

        %% Find the run onset times

        %Initialize variable
        run_start_times = [];

        %Loop over rows of log file
        for row = 1:length(DTA)

            %Search for string indicating start of a run
            if  any(~cellfun('isempty',strfind(DTA(row),'Keypress: 5'))) && any(~cellfun('isempty',strfind(DTA(row+1),'Exiting text screen'))) ...
                    && ~any(~cellfun('isempty',strfind(DTA(row+4),'Calibration'))) && ~any(~cellfun('isempty',strfind(DTA(row+5),'Starting trial 13')))...
                    && ~any(~cellfun('isempty',strfind(DTA(row+5),'Starting trial 14')))

                %Store the run start times
                run_start = DTA(row);
                run_start = run_start{1};
                run_start = run_start(1:7);
                run_start = str2num(run_start);

               %Create index of run start rows
               run_start_times = cat(1, run_start_times, run_start); 

            end

        end

        %% Find face, question, and button press events within trials

        %Within trial event detection - loop between rows of a trial
        for index = 1:length(trial_start_times)

            %Define the run number current trial belongs to
            if isequal(file, 1) %If running the first log file
                
                %Display trial number
                disp(['Running trial ', num2str(index)])
                
                %Note: for the first log file the run_number and
                %log_run_number are the same. For subsequent log files,
                %these two variables will diverge. 

                %Define run number by trial
                if index > 0 && index < 25

                    run_number = 1; %Absolute run number across log files
                    log_run_number = run_number; %Relative run number in current log file

                elseif index > 24 && index < 49

                    run_number = 2;
                    log_run_number = run_number;

                elseif index > 48 && index < 73

                    run_number = 3;
                    log_run_number = run_number;

                elseif index > 72 && index < 97

                    run_number = 4;
                    log_run_number = run_number;

                elseif index > 96 && index < 121

                    run_number = 5;
                    log_run_number = run_number;

                end

            elseif num_log_files > 1 %If running the 2+ log file

                %Define run number by trial relative to this log file
                if index > 0 && index < 25

                    log_run_number = 1;

                elseif index > 24 && index < 49

                    log_run_number = 2;

                elseif index > 48 && index < 73

                    log_run_number = 3;

                elseif index > 72 && index < 97

                    log_run_number = 4;

                elseif index > 96 && index < 121

                    log_run_number = 5;

                end

                %Add trials from previous log files to index to account for previous log files
                adjusted_index = index + num_trial_start_times;
                
                %Display trial number
                disp(['Running trial ', num2str(adjusted_index)])

                %Define run number by trial absolute for the session
                if adjusted_index > 0 && adjusted_index < 25

                    run_number = 1;

                elseif adjusted_index > 24 && adjusted_index < 49

                    run_number = 2;

                elseif adjusted_index > 48 && adjusted_index < 73

                    run_number = 3;

                elseif adjusted_index > 72 && adjusted_index < 97

                    run_number = 4;

                elseif adjusted_index > 96 && adjusted_index < 121

                    run_number = 5;

                end

            end

            %Make sure index does not exceed the total matrix size
            if index < length(trial_start_times)

                %Define last row of trial
                end_row = trial_start_times(index+1)-1;

            %If the last trial in a session
            elseif isequal(index, length(trial_start_times))

                %Define last row of trial
                end_row = length(DTA);

            end

            %Loop over rows of log file bounded within a trial
            for row = trial_start_times(index):end_row

                %Find blank faces and skip to next trial when found
                if any(~cellfun('isempty',strfind(DTA(row),'Drawing blank stimulus'))) && any(~cellfun('isempty',strfind(DTA(row+14),'Task Relevant'))) ...
                        || any(~cellfun('isempty',strfind(DTA(row),'Drawing blank stimulus'))) && any(~cellfun('isempty',strfind(DTA(row+15),'Task Relevant'))) ...
                        || any(~cellfun('isempty',strfind(DTA(row),'Drawing blank stimulus'))) && any(~cellfun('isempty',strfind(DTA(row+16),'Task Relevant')))

                    %Display as blank stimulus
                    disp('Blank sitmulus trial')

                    %End loop - Skip to next trial
                    break

                end

                %Find the time the face was drawn - if there is a face that is task relevant search for keypresses
                if any(~cellfun('isempty',strfind(DTA(row),'Drawing stimulus'))) && any(~cellfun('isempty',strfind(DTA(row+14),'Task Relevant'))) ...
                        || any(~cellfun('isempty',strfind(DTA(row),'Drawing stimulus'))) && any(~cellfun('isempty',strfind(DTA(row+15),'Task Relevant'))) ...
                        || any(~cellfun('isempty',strfind(DTA(row),'Drawing stimulus'))) && any(~cellfun('isempty',strfind(DTA(row+16),'Task Relevant')))

                    disp('Found face relevant stimulus')

                    %Save variable with face row number
                    face_row = row;

                    %Store the time the face appeared
                    current_face_time = DTA(face_row);
                    current_face_time = current_face_time{1};
                    current_face_time = current_face_time(1:7);
                    current_face_time = str2num(current_face_time);
                    all_relevant_face_times = cat(1,all_relevant_face_times, current_face_time); 

                    %Find the postion of the face
                    if any(~cellfun('isempty',strfind(DTA(face_row+3),'FACE: pos ='))) || any(~cellfun('isempty',strfind(DTA(face_row+4),'FACE: pos =')))

                       %Store the time the face appeared
                       if any(~cellfun('isempty',strfind(DTA(face_row+3),'FACE: pos =')))

                           %Define face position variable
                           current_relevant_face_pos = DTA(face_row+3);

                       elseif any(~cellfun('isempty',strfind(DTA(face_row+4),'FACE: pos =')))

                           %Define face position variable
                           current_relevant_face_pos = DTA(face_row+4);

                       end

                       %Store face positions in matrix
                       relevant_face_position = cat(1,relevant_face_position, current_relevant_face_pos);  

                       %Covert face position to location number
                       if contains(current_relevant_face_pos, '-0.8,  0.5') || contains(current_relevant_face_pos, ' 0. ,  0.4')

                           %Define as position 1
                           current_relevant_face_pos_num = 1;

                       elseif contains(current_relevant_face_pos, ' 0.8,  0.5') || contains(current_relevant_face_pos, ' 0. , -0.4')

                           %Define as position 2
                           current_relevant_face_pos_num = 2; 

                       elseif contains(current_relevant_face_pos, '-0.8, -0.5') || contains(current_relevant_face_pos, '-0.25,  0.')

                           %Define as position 3
                           current_relevant_face_pos_num = 3; 

                       elseif   contains(current_relevant_face_pos, ' 0.8, -0.5') || contains(current_relevant_face_pos, ' 0.25,  0.')

                           %Define as position 4
                           current_relevant_face_pos_num = 4; 

                       end

                    end

                    %Find the opacity of the face
                    if any(~cellfun('isempty',strfind(DTA(face_row+4),'FACE: opacity ='))) || any(~cellfun('isempty',strfind(DTA(face_row+5),'FACE: opacity =')))

                       %Store the time the face appeared
                       if any(~cellfun('isempty',strfind(DTA(face_row+4),'FACE: opacity =')))

                           %Define face position variable
                           current_relevant_face_opacity = DTA(face_row+4);

                       elseif any(~cellfun('isempty',strfind(DTA(face_row+5),'FACE: opacity =')))

                           %Define face position variable
                           current_relevant_face_opacity = DTA(face_row+5);

                       end

                       %Store face positions in matrix
                       relevant_face_opacity = cat(1,relevant_face_opacity, current_relevant_face_opacity); 

                       %Covert face position to location number
                       if contains(current_relevant_face_opacity, ' 1.0') 

                           %Define opaque
                           current_relevant_face_opacity = 'opaque';
                           disp(['Opacity ',current_relevant_face_opacity])

                       else

                           %Define threshold
                           current_relevant_face_opacity = 'threshold';
                           disp(['Opacity ',current_relevant_face_opacity])

                       end

                    end

                    %% Search for question and button press events post face presentation

                    %Loop over log file rows to find question periods
                    for row = face_row:end_row

                        % Find the perception question time
                        if any(~cellfun('isempty',strfind(DTA(row),'Asking perception question')))

                            %Display that the perception question was found
                            disp('Found perception question')

                            %Search for the keypress response after the perception question
                            for search = 1:100

                                %Make sure index does not exceed the total matrix size
                                if row+search > length(DTA)

                                    %Exit for loop
                                    break

                                end

                                %Find perception button press
                                if any(~cellfun('isempty',strfind(DTA(row+search),'Keypress: 1')))

                                   %Current perception keypress value
                                   current_per_keypress_num = 1;

                                   %Determine perception response 
                                   if isequal(Yes, 1)

                                      %Perceived stimulus
                                      perception = 1;

                                   elseif isequal(Yes, 2)

                                      %Did not perceive stimulus
                                      perception = 0;

                                   end

                                   %Store the time the face appeared
                                   current_per_keypress = DTA(row+search);
                                   current_per_keypress = current_per_keypress{1};
                                   current_per_keypress = current_per_keypress(1:7);
                                   current_per_keypress = str2num(current_per_keypress);
                                   all_per_ques_keypresses = cat(1,all_per_ques_keypresses, current_per_keypress);

                                   %Break out of loop 
                                   break

                                %Find perception button press
                                elseif any(~cellfun('isempty',strfind(DTA(row+search),'Keypress: 2')))

                                   %Current perception keypress value
                                   current_per_keypress_num = 2;

                                   %Determine perception response
                                   if isequal(Yes, 1)

                                      %Did not perceive stimulus
                                      perception = 0;

                                   elseif isequal(Yes, 2)

                                      %Perceived stimulus
                                      perception = 1;

                                   end

                                   %Store the time the face appeared
                                   current_per_keypress = DTA(row+search);
                                   current_per_keypress = current_per_keypress{1};
                                   current_per_keypress = current_per_keypress(1:7);
                                   current_per_keypress = str2num(current_per_keypress);
                                   all_per_ques_keypresses = cat(1,all_per_ques_keypresses , current_per_keypress);

                                   %Break out of loop 
                                   break

                                end

                            end

                        end      

                        %% Find the location question time
                        if any(~cellfun('isempty',strfind(DTA(row),'Asking location question')))

                            %Display that location question was found
                            disp('Found location question')

                            %Loop over rows
                            for search = 1:100

                                %Make sure index does not exceed the total matrix size
                                if row+search > length(DTA)

                                    %Break out of loop
                                    break

                                end

                                %Reset correct location logical
                                correct_location = 0; %Default is incorrect location

                                %Find perception button press
                                if any(~cellfun('isempty',strfind(DTA(row+search),'Keypress: 1')))

                                   %Current location keypress value
                                   current_loc_keypress_num = 1;

                                    %Check if location is correct
                                    if isequal(current_relevant_face_pos_num,1)

                                        %Set correct location logical to 1
                                        correct_location = 1;

                                    end

                                   %Store the time the face appeared
                                   current_loc_keypress = DTA(row+search);
                                   current_loc_keypress = current_loc_keypress{1};
                                   current_loc_keypress = current_loc_keypress(1:7);
                                   current_loc_keypress = str2num(current_loc_keypress);
                                   all_loc_ques_keypresses = cat(1,all_loc_ques_keypresses, current_loc_keypress);

                                   %Break out of loop 
                                   break

                                %Find perception button press
                                elseif any(~cellfun('isempty',strfind(DTA(row+search),'Keypress: 2')))

                                   %Current location keypress value
                                   current_loc_keypress_num = 2;

                                    %Check if location is correct
                                    if isequal(current_relevant_face_pos_num, 2)

                                        %Set correct location logical to 1
                                        correct_location = 1;

                                    end

                                   %Store the time the face appeared
                                   current_loc_keypress = DTA(row+search);
                                   current_loc_keypress = current_loc_keypress{1};
                                   current_loc_keypress = current_loc_keypress(1:7);
                                   current_loc_keypress = str2num(current_loc_keypress);
                                   all_loc_ques_keypresses = cat(1,all_loc_ques_keypresses, current_loc_keypress);

                                   %Break out of loop 
                                   break

                                %Find perception button press
                                elseif any(~cellfun('isempty',strfind(DTA(row+search),'Keypress: 3')))

                                   %Current location keypress value
                                   current_loc_keypress_num = 3;

                                    %Check if location is correct
                                    if isequal(current_relevant_face_pos_num,3)

                                        %Set correct location logical to 1
                                        correct_location = 1;

                                    end

                                   %Store the time the face appeared
                                   current_loc_keypress = DTA(row+search);
                                   current_loc_keypress = current_loc_keypress{1};
                                   current_loc_keypress = current_loc_keypress(1:7);
                                   current_loc_keypress = str2num(current_loc_keypress);
                                   all_loc_ques_keypresses = cat(1,all_loc_ques_keypresses, current_loc_keypress);

                                   %Break out of loop 
                                   break

                                %Find perception button press
                                elseif any(~cellfun('isempty',strfind(DTA(row+search),'Keypress: 4')))

                                   %Current location keypress value
                                   current_loc_keypress_num = 4;

                                    %Check if location is correct
                                    if isequal(current_relevant_face_pos_num,4)

                                        %Set correct location logical to 1
                                        correct_location = 1;

                                    end

                                   %Store the time the face appeared
                                   current_loc_keypress = DTA(row+search);
                                   current_loc_keypress = current_loc_keypress{1};
                                   current_loc_keypress = current_loc_keypress(1:7);
                                   current_loc_keypress = str2num(current_loc_keypress);
                                   all_loc_ques_keypresses = cat(1,all_loc_ques_keypresses, current_loc_keypress);

                                   %Break out of loop 
                                   break

                                end

                            end

                        end

                    end

                    %% Determine if CP or CnP Trial
                    
                    %Note: Use the relative run number (log_run_number) to
                    %specify the run_start_times and then add to the CP/CnP
                    %variable defined by the absolute run number. 
                    
                    %Perceived and location correct
                    if isequal(perception, 1) && isequal(correct_location, 1) && isequal(current_relevant_face_opacity, 'opaque')

                        %Save the CP face stimulus times - Subtract run start time
                        eval(['CP_opaque_face_times_run_',num2str(run_number),' = cat(1,CP_opaque_face_times_run_',num2str(run_number),', current_face_time-run_start_times(log_run_number));'])

                    elseif isequal(perception, 0) && isequal(correct_location, 0) && isequal(current_relevant_face_opacity, 'opaque') 

                        %Save the CnP face stimulus times - Subtract run start time
                        eval(['CnP_opaque_face_times_run_',num2str(run_number),' = cat(1,CnP_opaque_face_times_run_',num2str(run_number),', current_face_time-run_start_times(log_run_number));'])

                    elseif isequal(perception, 1) && isequal(correct_location, 1) && isequal(current_relevant_face_opacity, 'threshold')

                        %Save the CP face stimulus times - Subtract run start time
                        eval(['CP_thres_face_times_run_',num2str(run_number),' = cat(1,CP_thres_face_times_run_',num2str(run_number),', current_face_time-run_start_times(log_run_number));'])

                        %Not perceived and location incorrect
                    elseif isequal(perception, 0) && isequal(correct_location, 0) && isequal(current_relevant_face_opacity, 'threshold') 

                        %Save the CnP face stimulus times - Subtract run start time
                        eval(['CnP_thres_face_times_run_',num2str(run_number),' = cat(1,CnP_thres_face_times_run_',num2str(run_number),', current_face_time-run_start_times(log_run_number));'])

                        %Leave loop after face and other events have been identified
                        break

                    end

                end

            end

        end

        %Count the number of trials from previous runs - This variable is
        %used to adjust the index above when specifing the absolute run
        %number
        if isequal(file, 1)

            %Number of trial start times found in log file and previous log files
            num_trial_start_times = size(trial_start_times, 1);

        %If log file 2+    
        else 

            %Number of trial start times found in log file and previous log files
            num_trial_start_times = num_trial_start_times + size(trial_start_times, 1);

        end

    end
    
    %Calculate Sum of CP and CnP Trials
    CP_thres_num = length(CP_thres_face_times_run_1) + length(CP_thres_face_times_run_2) + length(CP_thres_face_times_run_3) + length(CP_thres_face_times_run_4) + length(CP_thres_face_times_run_5);
    CnP_thres_num = length(CnP_thres_face_times_run_1) + length(CnP_thres_face_times_run_2) + length(CnP_thres_face_times_run_3) + length(CnP_thres_face_times_run_4) + length(CnP_thres_face_times_run_5);
    CP_opaque_num = length(CP_opaque_face_times_run_1) + length(CP_opaque_face_times_run_2) + length(CP_opaque_face_times_run_3) + length(CP_opaque_face_times_run_4) + length(CP_opaque_face_times_run_5);
    CnP_opaque_num = length(CnP_opaque_face_times_run_1) + length(CnP_opaque_face_times_run_2) + length(CnP_opaque_face_times_run_3) + length(CnP_opaque_face_times_run_4) + length(CnP_opaque_face_times_run_5);

    %% Save data
    cd(savefolder)
    save MRI_CP_CnP_event_times.mat CP* CnP*

end