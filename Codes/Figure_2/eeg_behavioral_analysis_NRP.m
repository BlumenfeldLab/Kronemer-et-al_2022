%% EEG Behavioral Analysis Function - Visual No Report Paradigm 

%Behavioral analysis to identify perception, opacity, first stimulus, and
%location information for each face stimulus used to categorize and cut 
%trials later in preprocessing.

%NOTE: We should expect 120 trials if completed 5 runs, 48 trials in the
%first session of 2 runs and 72 trials in the second session of 3 runs.
%This is rarely different in some participants because of differences in
%how many runs were completed. 

%Written by: Sharif I. Kronemer
%Edited: 1/19/2021

function eeg_behavioral_analysis_NRP(Beh_dir, events_dir, relevant_location)

    %Create structure with the names of the raw data files across runs
    cd(Beh_dir)
    Beh_data = dir('*.csv'); 

    %% Run behavioral analysis to find and categorize trial types

    %Loop through the raw_files that equals the number of available behavioral run csv files
    for raw_files = 1:size(Beh_data,1)

        disp(['Running behavioral analysis - Session ', num2str(raw_files)])
        
        %Capture the import options
        cd(Beh_dir)
        options = detectImportOptions(Beh_data(raw_files,1).name);
        
        %Open Table - Utilize readtable function 
        table_raw_data = readtable(Beh_data(raw_files,1).name, options);
        
        %Capture the table column labels
        table_labels = table_raw_data.Properties.VariableNames;
        
        %Clear calibration trials and empty column variable
        try 
            
            if strcmp(table_raw_data.QUESTIONTYPE(1),'CALIBRATION')
            
                table_raw_data(strcmp(table_raw_data.QUESTIONTYPE,'CALIBRATION'),:) = [];
            
            end
        
        end
        
        %Remove last row of table if empty
        %try isempty(table_raw_data.QUESTIONTYPE{end});
        try 
            
            if isnan(table_raw_data.ABSOLUTETRIAL(end))
            
                %Delete
                table_raw_data(end,:) = [];
                
            end
            
        catch
            
            %Remove empty row
            if isempty(table_raw_data.TaskParadigm{end})

                %Delete
                table_raw_data(end,:) = [];

            end
            
        end

        %% Finds the data type of columns in raw file

        %Generate additional columns in table for new variables filled with NaN
        table_raw_data.PerceptionAccuracy = repmat({'NaN'}, size(table_raw_data,1),1);
        table_raw_data.LocationAccuracy = repmat({'NaN'}, size(table_raw_data,1),1);

        %% Finds the different perception types, reaction times, and stimulus hemifield

        %Loop through the rows of the raw file
        for row = 1:size(table_raw_data,1)
            
            if strcmp(relevant_location, 'Quadrant Relevant')
            
                % Fill in the perception accuracy for that row
                if strcmp(char(table_raw_data.QuadrantFaceShown(row)), 'True') && isequal(table_raw_data.PerceptionAnswer(row), 1) %Face present and perceived
                    table_raw_data.PerceptionAccuracy{row} = 'TP'; 

                elseif strcmp(char(table_raw_data.QuadrantFaceShown(row)), 'True') && isequal(table_raw_data.PerceptionAnswer(row), 0) %Face present and not perceived
                    table_raw_data.PerceptionAccuracy{row} = 'FN';
                    
                elseif strcmp(char(table_raw_data.QuadrantFaceShown(row)), 'False') && isequal(table_raw_data.PerceptionAnswer(row), 1) %Face absent and "percevied" 
                    table_raw_data.PerceptionAccuracy{row} = 'FP';
                    
                elseif strcmp(char(table_raw_data.QuadrantFaceShown(row)), 'False') && isequal(table_raw_data.PerceptionAnswer(row), 0) %Face absent and not perceived
                    table_raw_data.PerceptionAccuracy{row} = 'TN';
                    
                end

                % Fill in location accuracy for that row
                if strcmp(table_raw_data.PerceptionAccuracy(row),'TP') && strcmp(table_raw_data.QuadrantFaceLocation(row), table_raw_data.LocationAnswer(row)) == 1 %TP and location correct
                    table_raw_data.LocationAccuracy{row} = 'Confirmed perceived';

                elseif strcmp(table_raw_data.PerceptionAccuracy(row),'FN') && strcmp(table_raw_data.QuadrantFaceLocation(row), table_raw_data.LocationAnswer(row)) == 0 %FN and location incorrect
                    table_raw_data.LocationAccuracy{row} = 'Confirmed not perceived';

                elseif strcmp(table_raw_data.PerceptionAccuracy(row),'FN') && strcmp(table_raw_data.QuadrantFaceLocation(row), table_raw_data.LocationAnswer(row)) == 1 %FN and location correct
                    table_raw_data.LocationAccuracy{row} = 'Correct guess';

                elseif strcmp(table_raw_data.PerceptionAccuracy(row),'TP') && strcmp(table_raw_data.QuadrantFaceLocation(row), table_raw_data.LocationAnswer(row)) == 0 %TP and location incorrect
                    table_raw_data.LocationAccuracy{row} = 'False perception';

                elseif strcmp(table_raw_data.PerceptionAccuracy(row),'TN') 
                    table_raw_data.LocationAccuracy{row} = 'True Negative';

                elseif strcmp(table_raw_data.PerceptionAccuracy(row),'FP') 
                    table_raw_data.LocationAccuracy{row} = 'False Positive';

                end 
                
            elseif strcmp(relevant_location, 'Center Relevant')
                
                % Fill in the perception accuracy for that row
                if strcmp(char(table_raw_data.CenterFaceShown(row)), 'True') && isequal(table_raw_data.PerceptionAnswer(row), 1) %Face present and perceived
                    table_raw_data.PerceptionAccuracy{row} = 'TP'; 

                elseif strcmp(char(table_raw_data.CenterFaceShown(row)), 'True') && isequal(table_raw_data.PerceptionAnswer(row), 0) %Face present and not perceived
                    table_raw_data.PerceptionAccuracy{row} = 'FN';
                    
                elseif strcmp(char(table_raw_data.CenterFaceShown(row)), 'False') && isequal(table_raw_data.PerceptionAnswer(row), 1) %Face absent and "percevied" 
                    table_raw_data.PerceptionAccuracy{row} = 'FP';
                    
                elseif strcmp(char(table_raw_data.CenterFaceShown(row)), 'False') && isequal(table_raw_data.PerceptionAnswer(row), 0) %Face absent and not perceived
                    table_raw_data.PerceptionAccuracy{row} = 'TN';
                    
                end

                % Fill in location accuracy for that row
                if strcmp(table_raw_data.PerceptionAccuracy(row),'TP') && strcmp(table_raw_data.CenterFaceLocation(row), table_raw_data.LocationAnswer(row)) == 1 %TP and location correct
                    table_raw_data.LocationAccuracy{row} = 'Confirmed perceived';

                elseif strcmp(table_raw_data.PerceptionAccuracy(row),'FN') && strcmp(table_raw_data.CenterFaceLocation(row), table_raw_data.LocationAnswer(row)) == 0 %FN and location incorrect
                    table_raw_data.LocationAccuracy{row} = 'Confirmed not perceived';

                elseif strcmp(table_raw_data.PerceptionAccuracy(row),'FN') && strcmp(table_raw_data.CenterFaceLocation(row), table_raw_data.LocationAnswer(row)) == 1 %FN and location correct
                    table_raw_data.LocationAccuracy{row} = 'Correct guess';

                elseif strcmp(table_raw_data.PerceptionAccuracy(row),'TP') && strcmp(table_raw_data.CenterFaceLocation(row), table_raw_data.LocationAnswer(row)) == 0 %TP and location incorrect
                    table_raw_data.LocationAccuracy{row} = 'False perception';

                elseif strcmp(table_raw_data.PerceptionAccuracy(row),'TN') 
                    table_raw_data.LocationAccuracy{row} = 'True Negative';

                elseif strcmp(table_raw_data.PerceptionAccuracy(row),'FP') 
                    table_raw_data.LocationAccuracy{row} = 'False Positive';

                end 
                
            end
       
        end 

        %% Generate Trial Identifiers
        
        %CP, CnP, Correct guess, etc.
        perception_identifier = table_raw_data.LocationAccuracy(:);
        
        %Center Relevant
        if strcmp(relevant_location, 'Center Relevant')
        
            %Location of face on screen
            relevant_face_location = table_raw_data.CenterFaceLocation(:); 
            irrelevant_face_location = table_raw_data.QuadrantFaceLocation(:);
            
            %See if column is a cell array (which will require a
            %transformation to double)
            if iscell(table_raw_data.CenterFaceOpacity(:))
                
                %Relevant face opacity
                relevant_face_opacity = str2double(table_raw_data.CenterFaceOpacity(:));
                
                %Irrelevant face opacity
                irrelevant_face_opacity = str2double(table_raw_data.QuadrantFaceOpacity(:));
            
            else
                
                %Relevant face opacity
                relevant_face_opacity = table_raw_data.CenterFaceOpacity(:);
                
                %Irrelevant face opacity
                irrelevant_face_opacity = table_raw_data.QuadrantFaceOpacity(:);
            
            end
        
        %Quadrant Relevant
        elseif strcmp(relevant_location, 'Quadrant Relevant')
            
            %Location of face on screen
            relevant_face_location = table_raw_data.QuadrantFaceLocation(:);
            irrelevant_face_location = table_raw_data.CenterFaceLocation(:);

            %See if column is a cell array (which will require a
            %transformation to double)
            if iscell(table_raw_data.CenterFaceOpacity(:))
            
                %Relevant face opacity
                relevant_face_opacity = str2double(table_raw_data.QuadrantFaceOpacity(:));

                %Irrelevant face opacity
                irrelevant_face_opacity = str2double(table_raw_data.CenterFaceOpacity(:));

            else
                
                %Relevant face opacity
                relevant_face_opacity = table_raw_data.QuadrantFaceOpacity(:);

                %Irrelevant face opacity
                irrelevant_face_opacity = table_raw_data.CenterFaceOpacity(:);
                
            end
            
        end
        
        %Find the first stimulus per trial and store information
        first_stim_identifier = table_raw_data.FirstStimulus(:); 

        %Save identifiers 
        cd(events_dir)
        save(['Raw_file_',num2str(raw_files),'_trial_identifiers.mat'], 'perception_identifier', 'relevant_face_location', 'irrelevant_face_location', 'first_stim_identifier', 'relevant_face_opacity', 'irrelevant_face_opacity')   
   
    end
    
end