%Identify Task Irrelevant Stimuli from Behavioral Log File

%The goal of this code is to load the behavioral log file and identify CP
%and CnP trials from the events registered with the log file. This code was
%particularly designed for the visual no report paradigm version 1.12. 

%Note: This code has been modified to allow managing multiple log files
%that are created when the behavioral task is restarted during a testing
%session.

%Written by: Sharif I. Kronemer
%Date: 7/23/2019
%Modified: 1/5/2021

function Find_Task_Irrelevant_Events_From_Log_File_Function(beh_folder, savefolder)

    %Load the behavioral log file
    cd(beh_folder)
    
    %Find the log files in behavioral directory
    log_files = dir('*.log'); 
    
    %Count the number of log files
    num_log_files = size(log_files,1); 

    %% Initialize variables

    all_opaque_irrelevant_face_times = [];
    all_threshold_irrelevant_face_times = [];
    all_blank_irrelevant_face_times = [];

    %Opaque run variables
    opaque_face_times_run_1 = [];
    opaque_face_times_run_2 = [];
    opaque_face_times_run_3 = [];
    opaque_face_times_run_4 = [];
    opaque_face_times_run_5 = [];

    %Threshold run variables
    threshold_face_times_run_1 = [];
    threshold_face_times_run_2 = [];
    threshold_face_times_run_3 = [];
    threshold_face_times_run_4 = [];
    threshold_face_times_run_5 = [];

    %Blank run variables
    blank_face_times_run_1 = [];
    blank_face_times_run_2 = [];
    blank_face_times_run_3 = [];
    blank_face_times_run_4 = [];
    blank_face_times_run_5 = [];
       
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

            %Search for string indicating start of a trial
            if  any(~cellfun('isempty',strfind(DTA(row),'Keypress: 5'))) && any(~cellfun('isempty',strfind(DTA(row+1),'Exiting text screen'))) ...
                    && ~any(~cellfun('isempty',strfind(DTA(row+4),'Calibration'))) && ~any(~cellfun('isempty',strfind(DTA(row+5),'Starting trial 13'))) ...
                    && ~any(~cellfun('isempty',strfind(DTA(row+5),'Starting trial 14')))

                %Store the run start times
                run_start = DTA(row);
                run_start = run_start{1};
                run_start = run_start(1:7);
                run_start = str2num(run_start);

               %Create index of trial start rows
               run_start_times = cat(1, run_start_times, run_start); 

            end

        end

        %% Find irrelevant face events within trials

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

            %Count number of blank stim
            blank_stim = 0;

            %Loop over rows of log file bounded within a trial
            for row = trial_start_times(index):end_row

                %Find blank faces
                if  any(~cellfun('isempty',strfind(DTA(row),'Drawing blank stimulus'))) && any(~cellfun('isempty',strfind(DTA(row+14),'Task Irrelevant'))) ...
                        || any(~cellfun('isempty',strfind(DTA(row),'Drawing blank stimulus'))) && any(~cellfun('isempty',strfind(DTA(row+15),'Task Irrelevant'))) ...
                        || any(~cellfun('isempty',strfind(DTA(row),'Drawing blank stimulus'))) && any(~cellfun('isempty',strfind(DTA(row+16),'Task Irrelevant')))

                    %Display as blank stimulus
                    disp('Blank sitmulus trial')
                    blank_stim = blank_stim +1;

                    %Save variable with face row number
                    face_row = row;

                    %Store the time the face appeared
                    current_face_time = DTA(face_row);
                    current_face_time = current_face_time{1};
                    current_face_time = current_face_time(1:7);
                    current_face_time = str2num(current_face_time);

                    %Store all opaque face times in one matrix 
                    all_blank_irrelevant_face_times = cat(1,all_blank_irrelevant_face_times, current_face_time); 

                    %Note: Use the relative run number (log_run_number) to
                    %specify the run_start_times and then add to the CP/CnP
                    %variable defined by the absolute run number. 
                    
                    %Save face time in matrix; Subtract current face tie with
                    %run start time
                    eval(['blank_face_times_run_',num2str(run_number),' = cat(1,blank_face_times_run_',num2str(run_number),', current_face_time-run_start_times(log_run_number));'])

                end

                %Find the time the face was drawn - if there is a face that is task relevant search for keypresses
                if any(~cellfun('isempty',strfind(DTA(row),'Drawing stimulus'))) && any(~cellfun('isempty',strfind(DTA(row+14),'Task Irrelevant'))) ...
                        || any(~cellfun('isempty',strfind(DTA(row),'Drawing stimulus'))) && any(~cellfun('isempty',strfind(DTA(row+15),'Task Irrelevant'))) ...
                        || any(~cellfun('isempty',strfind(DTA(row),'Drawing stimulus'))) && any(~cellfun('isempty',strfind(DTA(row+16),'Task Irrelevant')))

                    disp('Found face irrelevant stimulus')

                    %Save variable with face row number
                    face_row = row;

                    %Store the time the face appeared
                    current_face_time = DTA(face_row);
                    current_face_time = current_face_time{1};
                    current_face_time = current_face_time(1:7);
                    current_face_time = str2num(current_face_time);          

                    %Determine the opacity of the face
                    if any(~cellfun('isempty',strfind(DTA(row+4),'opacity = 1.0')))||any(~cellfun('isempty',strfind(DTA(row+5),'opacity = 1.0')))

                        disp('Fully opaque irrelevant face')

                        %Store all opaque face times in one matrix 
                        all_opaque_irrelevant_face_times = cat(1,all_opaque_irrelevant_face_times, current_face_time); 

                        %Note: Use the relative run number (log_run_number) to
                        %specify the run_start_times and then add to the CP/CnP
                        %variable defined by the absolute run number. 
                        
                        %Save face time in matrix; Subtract current face tie with
                        %run start time
                        eval(['opaque_face_times_run_',num2str(run_number),' = cat(1,opaque_face_times_run_',num2str(run_number),', current_face_time-run_start_times(log_run_number));'])
                    
                    %Determine the opacity of the face - special case for
                    %569 Center Relevant where the quadrant opacity in run
                    %1 was above 1.
                    elseif any(~cellfun('isempty',strfind(DTA(row+4),'opacity = 13')))||any(~cellfun('isempty',strfind(DTA(row+5),'opacity = 13')))

                        disp('Fully opaque irrelevant face')
                        
                        %Store all opaque face times in one matrix 
                        all_opaque_irrelevant_face_times = cat(1,all_opaque_irrelevant_face_times, current_face_time); 

                        %Note: Use the relative run number (log_run_number) to
                        %specify the run_start_times and then add to the CP/CnP
                        %variable defined by the absolute run number. 
                        
                        %Save face time in matrix; Subtract current face tie with
                        %run start time
                        eval(['opaque_face_times_run_',num2str(run_number),' = cat(1,opaque_face_times_run_',num2str(run_number),', current_face_time-run_start_times(log_run_number));'])
                                    
                    else

                        disp('Threshold irrelevant face')

                        %Store all threshold face times in one matrix
                        all_threshold_irrelevant_face_times = cat(1,all_threshold_irrelevant_face_times, current_face_time); 

                        %Note: Use the relative run number (log_run_number) to
                        %specify the run_start_times and then add to the CP/CnP
                        %variable defined by the absolute run number. 
                    
                        %Save face time in matrix; Subtract current face time with
                        %run start time
                        eval(['threshold_face_times_run_',num2str(run_number),' = cat(1,threshold_face_times_run_',num2str(run_number),', current_face_time-run_start_times(log_run_number));'])

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

    %Calculate Sum of Irrelevant Trials
    blank_face_num = length(blank_face_times_run_1) + length(blank_face_times_run_2) + length(blank_face_times_run_3) + length(blank_face_times_run_4) + length(blank_face_times_run_5);
    opaque_face_num = length(opaque_face_times_run_1) + length(opaque_face_times_run_2) + length(opaque_face_times_run_3) + length(opaque_face_times_run_4) + length(opaque_face_times_run_5);
    threshold_face_num = length(threshold_face_times_run_1) + length(threshold_face_times_run_2) + length(threshold_face_times_run_3) + length(threshold_face_times_run_4) + length(threshold_face_times_run_5);

    %% Save data
    cd(savefolder)
    save MRI_irrelevant_face_event_times.mat threshold* opaque* blank_face* all_opaque_irrelevant_face_times all_threshold_irrelevant_face_times all_blank_irrelevant_face_times run_start_times trial_start_times

end