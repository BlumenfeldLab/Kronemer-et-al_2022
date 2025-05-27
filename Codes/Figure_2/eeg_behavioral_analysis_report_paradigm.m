%% EEG Behavioral Analysis Code for Visual Perception Task

%***Inputs***
%Subject ID
%Behavioral data directory
%Save data directory

%***Outputs***
%Events directory where the behavioral analysis is saved
%Behavioral data directories 

%Written by: Sharif I. Kronemer
%Edited: 4/29/2021
%Edited: Julia Ding 11/15/2019

function eeg_behavioral_analysis_report_paradigm(ID, Beh_dir, save_dir, events_dir)

    %Navigate to behavioral directory
    cd(Beh_dir)

    %Create a variable with the names of folders inside the Beh_dir directory 
    Beh_folders = struct2table(dir(Beh_dir));
    
    %Initialize folder variable
    folders = [];
    
    %Loop over folders
    for i = 1:size(Beh_folders,1)
        
        %Get name of folder
        name = Beh_folders.name(i);
        
        %Find folders with Run in name
        if contains(name,'Run') & ~contains(name, 'caliFile')
            folders = [folders; name];
            
        end 
        
    end
    
    %Find the number of behavioral folders
    num_beh_files = length(folders);
    
    %Individual subject correction
    if isequal(ID, '243EG')
        
        num_beh_files = 2;
        
    end
    
    %Create a cell with all the names of the csv file for behavioral data
    Beh_data = {};

    %For all subjects
    for n = 1:num_beh_files 

        cd(string(fullfile(Beh_dir, folders(n)))); %Enter run folder with the EEG data
        Beh_data{n} = dir('*.csv'); %Create structure with the names of the raw data files across runs

    end

    %Set Exception for 336MA that has all behavioral files but missing 1 raw file
    if isequal(ID, '366MA')

       Beh_data = Beh_data(1,2); %Set the Beh_data cell to only the 3rd csv file with runs 3 and 4
       folders(1) = []; %Remove Calibration row from folder directory so that raw files pulls out the correct file for this person

    end

    %% Create and define save file directory 

%     %Check for the existences of the epochs directory
%     if exist(fullfile(save_dir, 'Epochs')) == 0
% 
%         %Create folder
%         mkdir(fullfile(save_dir, 'Epochs')); %Creates directory for epochs if it does not already exist
% 
%     end     

%     %Create a directory for saving Psychopy and DIN Marker event times if it does not exist
%     if exist(fullfile(save_dir, 'Event Times and Identifiers')) == 0
% 
%         %Create folder
%         mkdir(fullfile(save_dir, 'Event Times and Identifiers'));
% 
%     end
% 
%     %Create a directory for saving Psychopy and DIN Marker event times if it does not exist
%     if exist(fullfile(save_dir, 'Sessions')) == 0
% 
%         %Create folder
%         mkdir(fullfile(save_dir, 'Sessions'));
% 
%     end

    %% Run behavioral analysis to find and categorize trial types

    %Loop through the raw_files that equals the number of available behavioral run csv files
    for raw_files = 1:size(Beh_data,2)

        disp(['Running behavioral analysis for session ', num2str(raw_files)])

        %Enter the .csv files into Behavior_data cell
        cd(string(fullfile(Beh_dir, folders(raw_files))))%Enter folder with appropriate csv file

        %% Find PsychoPy Times CP, CnP, Button Press, Trial onset, and Question times

        %Utilize readtable function    
        current_structure = Beh_data{1,raw_files};
        name_of_excel_file = current_structure.name; 
        table_raw_data = readtable(name_of_excel_file);
        table_labels = {'BLOCK NUMBER' 'Trial start time' 'TRIAL TYPE' 'ABSOLUTE TRIAL' 'QUESTION TYPE' 'TRIAL' 'Delay' 'Jitter prestimulus delay value' 'FaceDrawStart' 'Face duration' 'Face opacity' 'Face quadrant' 'Face shown' 'Actual prestimulus delay (face presentation time minus trial start)' 'Perception question time'	'Perception keypress' 'Perception keypress time' 'Perception answer' 'Location question time' 'Location keypress' 'Location keypress time' 'Location answer' 'Trial Tag'};

        %Cut extra column
        if isequal(size(table_raw_data,2), 24)
            
           %Remove last column
           table_raw_data = table_raw_data(:,1:end-1);

        end

        %Define raw variable - add labels to columns
        raw = [table_labels; table2cell(table_raw_data)];
        
        %Replace true or false strings with binary values
        if sum(ismember(raw(:,13), 'True')) > 0

            tru_idx = ismember(raw(:,13), 'True');
            fal_idx = ismember(raw(:,13), 'False');

        elseif  sum(ismember(raw(:,13), 'TRUE')) > 0

            tru_idx = ismember(raw(:,13), 'TRUE');
            fal_idx = ismember(raw(:,13), 'FALSE');

        end
        
        %Replace true/false with 1/0
        raw(tru_idx,13) = num2cell(1);
        raw(fal_idx, 13) = num2cell(0);

        %Remove data associated with calibration trials 
        calibration_row = [];

        %Finds column indicating trial type
        trialtype = find(strcmp('TRIAL TYPE',raw(1,:))); 

        %Loop through the number of rows available in raw
        for i = 1:size(raw,1)

            if strfind(raw{i,trialtype},'CALIBRATION') %Finds calibration data
                
                %Create a matrix of calibration rows
                calibration_row = [calibration_row i];

            end

        end

        %Empty all the rows that are related to calibration in raw
        raw(calibration_row,:) = [];

        % Skip if the behavioral file is only calibration 
        if size(raw,1) == 1

           disp('Skipping analysis of a calibration behavioral file.')

        end

        %If the first row is blank, delete it.
        if isnan(raw{2,1})
        
            raw(2,:) = []; 
        
        end

        %Find number of rows and columns
        row = size(raw,1); 
        column = size(raw,2); 

        %% Finds the data type of columns in raw file

        %Generate additional columns in raw for new variables
        raw{1,column+1} = 'Perception Accuracy';
        raw{1,column+2} = 'Perception RT';
        raw{1,column+3} = 'Location Accuracy';
        raw{1,column+4} = 'Location RT';
        raw{1,column+5} = 'Hemifield';

        %Loop through the number of columns available in the raw file
        for i = 1:size(raw,2)

            if strcmp('TRIAL TYPE',raw{1,i})
                trialtype = i;
            elseif strcmp('Delay',raw{1,i})
                delay = i;
            elseif strcmp('Perception answer',raw{1,i})
                perceptionanswer = i;
            elseif strcmp('Perception Accuracy',raw{1,i})
                perception = i;
            elseif strcmp('Perception question time',raw{1,i})
                questionp = i;       
            elseif strcmp('Location answer',raw{1,i})
                locationanswer = i;
            elseif strcmp('Location Accuracy',raw{1,i})
                location = i;
            elseif strcmp('Location question time',raw{1,i})
                questionl = i; 
            elseif strcmp('FaceDrawStart',raw{1,i})
                facedraw = i;
            elseif strcmp('BLOCK NUMBER',raw{1,i})
                block = i;
            elseif strcmp('Perception keypress time',raw{1,i})
                buttonp = i;
            elseif strcmp('Location keypress time',raw{1,i})
                buttonl = i;
            elseif strcmp('Perception RT',raw{1,i})
                reactiontimep = i;
            elseif strcmp('Location RT',raw{1,i})
                reactiontimel = i;
            elseif strcmp('Trial start time',raw{1,i})
                trialstart = i;
            elseif strcmp('TRIAL',raw{1,i})
                trial = i;
            elseif strcmp('Face shown',raw{1,i})
                faceshown = i;
            elseif strcmp('Face quadrant',raw{1,i})
                quadrant = i;    
            elseif strcmp('Face opacity',raw{1,i})
                opacity = i;
            elseif strcmp('Hemifield',raw{1,i})
                hemifield = i;       
            end

        end

        %Finds the total number of blocks completed (finds the max value in Block column of raw)
        totalblocks = max(cell2mat(raw(2:end,block)));

        %% Finds the different perception types, reaction times, and stimulus hemifield

        %Loop through the rows of the raw file
        for j = 2:size(raw,1) %Starts from row 2 in order to omit the label row

            %Use only the rows with actual run trials, omit calibration rows - by finding only trials that have NOISE or MOVIE in the trial type column 
            if strcmp(raw{j,trialtype},'NOISE') || strcmp(raw{j,trialtype},'MOVIE')

                % Fill in the perception accuracy for that row (j)
                if raw{j,faceshown} == 1 && raw{j,perceptionanswer} == 1 %Face present and perceived
                    raw{j,perception} = 'TP'; %True Positive

                elseif raw{j,faceshown} == 1 && raw{j,perceptionanswer} == 0 %Face present and not perceived
                    raw{j,perception} = 'FN'; %False Negative

                elseif raw{j,faceshown} == 0 && raw{j,perceptionanswer} == 1 %Face absent and "percevied" 
                    raw{j,perception} = 'FP'; %False Positive

                elseif raw{j,faceshown} == 0 && raw{j,perceptionanswer} == 0 %Face absent and not perceived
                    raw{j,perception} = 'TN'; %True Negative

                end

                % Fill in perception reaction time for that row (j)
                raw{j,reactiontimep} = raw{j,buttonp} - raw{j,questionp};

                % Fill in location accuracy for that row (j)
                if strcmp(raw{j,perception},'TP') && strcmp(raw{j,quadrant},raw{j,locationanswer}) == 1 %TP and location correct
                    raw{j,location} = 'Confirmed perceived';

                elseif strcmp(raw{j,perception},'FN') && strcmp(raw{j,quadrant},raw{j,locationanswer}) == 0 %FN and location incorrect
                    raw{j,location} = 'Confirmed not perceived';

                elseif strcmp(raw{j,perception},'FN') && strcmp(raw{j,quadrant},raw{j,locationanswer}) == 1 %FN and location correct
                    raw{j,location} = 'Correct guess';

                elseif strcmp(raw{j,perception},'TP') && strcmp(raw{j,quadrant},raw{j,locationanswer}) == 0 %TP and location incorrect
                    raw{j,location} = 'False perception';

                elseif strcmp(raw{j,perception},'TN') 
                    raw{j,location} = 'True Negative';

                elseif strcmp(raw{j,perception},'FP') 
                    raw{j,location} = 'False Positive';

                end              

                % Fill in location reaction time for that row (j)
                raw{j,reactiontimel} = raw{j,buttonl} - raw{j,questionl}; %Onset of location question minus button press time

                % Fill in the left vs right hemifield of the stimulus
                if strcmp(raw{j,quadrant},'[-0.8, 0.5]')
                    raw{j,hemifield} = 'Top Left';

                elseif strcmp(raw{j,quadrant},'[-0.8, -0.5]') 
                    raw{j,hemifield} = 'Bottom Left';

                elseif strcmp(raw{j,quadrant},'[0.8, 0.5]') 
                    raw{j,hemifield} = 'Top Right';

                elseif strcmp(raw{j,quadrant},'[0.8, -0.5]')
                    raw{j,hemifield} = 'Bottom Right';

                end

            end
        end 

        %% Generate a list of times when faces, question, buttons appear in Psychopy    

        %Create a vector of NaNs that is the length of the number of events of that type
        psychopy_faces = NaN(1,size(raw,1)-1); %1 column x (rows in raw -1)
        psychopy_questions = NaN(2,size(raw,1)-1); %2 column x (rows in raw -1)
        psychopy_buttons = NaN(2,size(raw,1)-1); %2 column x (rows in raw -1)

        %Populates the Psychopy variables (e.g., psychopy_faces, etc.) with the times of events
        for w = 1:size(raw,1)-1

            psychopy_faces(1,w) = raw{w+1,facedraw}; %Add time of face onset from the face draw column

            psychopy_questions(1,w) = raw{w+1,questionp}; %Add time of the perception question
            psychopy_questions(2,w) = raw{w+1,questionl}; %Add time of the location question

            psychopy_buttons(1,w) = raw{w+1,buttonp}; %Add time of the perception button press
            psychopy_buttons(2,w) = raw{w+1,buttonl}; %Add time of teh location button press

        end

        %% Generate Trial Identifiers (e.g., CP, CnP, correct guess, etc.)
        
        %CP, CnP, Correct guess, etc.
        perception_identifier = raw(2:end,location); 

        %Location of face on screen
        hemifield_identifier = raw(2:end,hemifield); 

        %Background of screen (movie or noise)
        background_identifier = raw(2:end,trialtype);

        %Post-stimulus duration (1 or 15s)
        delay_identifier = cell2mat(raw(2:end,delay)); 

        %Save identifiers 
        cd(events_dir)
        save(['Raw_file_',num2str(raw_files+1),'_trial_identifiers.mat'], 'perception_identifier', 'hemifield_identifier', 'background_identifier', 'delay_identifier')   

    end
    
end