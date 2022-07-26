%% Find Face Pulse Times in RNS Data and Cut RNS Epochs

%This code loads RNS sessions, find marks, and cuts epochs

% Written by: Sharif I. Kronemer
% Date: 9/14/2019
% Modified: 5/17/2021

clear

%% Prompts

prompt = 'Subject ID: '; 
ID = input(prompt, 's');

%% Directories and Paths

%Paths
addpath('R:\RNS Study\Analysis\Analysis Code\MATLAB function for reading DAT files')
addpath('R:\RNS Study\Analysis\Analysis Code')

%Directories
data_dir = ['R:\RNS Study\Subject Raw Data\',ID,'\RNS Data'];
beh_dir = ['R:\RNS Study\Subject Raw Data\',ID,'\Behavioral Data\Calibration and All Runs'];
%beh_results = ['P:\RNS Study\Analysis\Subject Analysis\',ID,'\Scalp EEG Analysis'];
beh_results = ['R:\RNS Study\Analysis\Subject Analysis\',ID,'\Scalp EEG Analysis\Event Times and Identifiers'];
rns_results = ['R:\RNS Study\Analysis\Subject Analysis\',ID,'\RNS Analysis'];

%Epochs dir
epoch_results = fullfile(rns_results,'Voltage Epochs');

%Make save folder if necessary
if ~exist(epoch_results)

   %Make folder if missing 
   mkdir(epoch_results)

end   

%Find run folders
cd(data_dir)
run_folders = dir('Run*');

%Initialize variable
run_list ={};

%Loop over directory structure
for row = 1:length(run_folders)

    %Find folder names
    run_list{row,1} = run_folders(row).name;

end

%Organize folder names in correct order
run_list = sort_nat(run_list);

%% Parameters

%Sampling rate
sampling_rate = 250;

%Epoch duration
epoch_duration = 2;

%RNS marker and trial onset delay
RNS_delay = 0.51;

%% Behavioral Data 

%Find csv files
cd(beh_dir)
csv_file = dir('*.csv');

% Load the relevant file
disp(['Load behavioral data'])
[num,text,raw] = xlsread(fullfile(beh_dir,csv_file.name));

% Remove calibration rows
raw(2:402,:) = [];
   
%Find number of rows and columns
row = length(raw(:,1));
column = length(raw(1,:));

 %% Find the appropriate columns and save to variables
 
%Loop over rows of the raw file
for i = 1:size(raw,2)

    %Find columns of interest
    if strcmp('TRIAL TYPE',raw{1,i})
        trialtype = i;
    
    elseif strcmp('Actual prestimulus delay (face presentation time minus trial start)',raw{1,i})
        prestim = i;
    
    elseif strcmp('BLOCK NUMBER',raw{1,i})
        block_num = i;
        
    elseif strcmp('ABSOLUTE TRIAL',raw{1,i})
        absolute_trial = i;
        
    end
    
end

%% Extract Event Times

%Initialize variable
all_CP_epochs = [];
all_CnP_epochs = [];

% Loop over runs
for run = 1:length(run_list)
    
    %% Behavioral Analysis
    
    %Load results
    load(fullfile(beh_results,'ttl_times_and_identity.mat'))
    
    %Find prestim times (+1 to ignore column text row)
    run_pre_stim = cell2mat(raw(find(cell2mat(raw(2:end,block_num)) == run)+1,prestim));
    
    %Trial identifier
    trial_ID = trial_identifiers(find(cell2mat(raw(2:end,block_num)) == run));
    
    %Select current run
    run_num = cell2mat(run_list(run));
    
    %Create run directory
    run_data_dir = fullfile(data_dir, run_num);
    
    %Find .dat files in run folder
    dat_files = dir(run_data_dir);
    dat_files(1:2) = [];
    
    %Loop over files
    for file = 1:length(dat_files)
        
        disp(['Loading data file ', num2str(file)])
        
        %Load data - added 0.765 multiplier to ReadECoGData function along
        %with 512 subtraction
        [ECoG_hdr, ECoG_data] = ReadECoGData(fullfile(run_data_dir,dat_files(file).name), 'test.lay');
        
        %Rename files
        eval(['ECoG_data_',num2str(file),' = ECoG_data;'])
        eval(['ECoG_hdr_',num2str(file),' = ECoG_data;'])
        
    end
    
    %Combine dat files 
    if isequal(ID, '471MH')
          
        %Special baselining procedure Note: channels 1,2 are subtracted 5 and channels 3,4 are subtracted 5.5    
        ECoG_data{1} = ECoG_data{1}-5;
        ECoG_data{2} = ECoG_data{2}-5;
        ECoG_data{3} = ECoG_data{3}-5.5;
        ECoG_data{4} = ECoG_data{4}-5.5;
        
        %No files to combine - Note: 300M only has 4 minutes of live ECoG
        %so the runs were truncated to 4 minutes.
        ECoG_data_combined = ECoG_data;
    
    else
        
        %All other subjects have files to be combined within runs - Note:
        %320M only has 8 minutes of live ECoG so the runs were truncated
        %to 8 minutes but the two continguous 4 minutes files are generated 
        ECoG_data_combined = {cat(2,ECoG_data_1{1,1},ECoG_data_2{1,1}), cat(2,ECoG_data_1{1,2},ECoG_data_2{1,2}),...
            cat(2,ECoG_data_1{1,3},ECoG_data_2{1,3}), cat(2,ECoG_data_1{1,4},ECoG_data_2{1,4})};
    
    end
    
    %% Find RNS Marks
    
    %Z-score data
    zscore_data = zscore(ECoG_data_combined{1,1});
    
    %Find outlier events corresponding with marks
    telemetry_times = find(zscore_data<-5)';
    
    %Individual subject corrections
    if isequal(ID, '491GS') && isequal(run, 3)
       
       %Delete markers during stimulation time of Run 3
       telemetry_times(527:541) = [];
       
    elseif isequal(ID, '536JD') && isequal(run, 8)
        
       %Delete markers from dropped signal
       telemetry_times(373:1586) = [];
       
    end
    
    %Update marks for 535BP because of dropped marks
    if isequal(ID, '535BP') || isequal(ID, '536JD')
       
        %Find difference between telemetry times
        telemetry_diff = diff(telemetry_times);
        
        %Loop for a 1,3,1,7,1 sequence in telemetry diff
        sequence_1 = strfind(telemetry_diff', [1,3,1,7,1]); 
        
        %Loop for a 1,7,1 sequence in telemetry diff
        sequence_2 = strfind(telemetry_diff',[1,7,1])-2; %Subtract 2, to find the onset of the 1,3,1,7,1 sequence
        
        %Find corresponding indices for sequence 1 and 2; the goal to
        %isolate the sequences that neglect the initial 1,3
        sequence_index = ~ismember(sequence_2,sequence_1)';

        %Add two values back subtracted before (to handle negative index
        %values if the first mark is affected in a trial)
        sequence_2 = sequence_2+2; 

        %Add the missing telemetry times back
        new_telemetry_times = telemetry_times; 
        new_telemetry_times = sort([telemetry_times; telemetry_times(sequence_2(sequence_index))-4; telemetry_times(sequence_2(sequence_index))-3;]);
        
        %Diff of new telemetry times
        new_tele_diff = diff(new_telemetry_times);
        
    end
    
    %Plot mark events with data
    
    %Loop over channels 
    for channel = 1:4
        
        RNS_figure = figure; 
        hold on
        title([ID, ' - Run ', num2str(run), ' - Channel ', num2str(channel)])
        ylabel('Voltage (uV)')
        xlabel('Samples')
        ylim([-500 100])

        %Plot figure
        plot(ECoG_data_combined{1,channel},'Marker','*','MarkerIndices',final_telemetry_times,'MarkerEdgeColor','r'); 

        %Make save folder if necessary
        if ~exist(fullfile(rns_results, 'Telemetry RNS Figures'))

           %Make folder if missing 
           mkdir(fullfile(rns_results, 'Telemetry RNS Figures'))

        end   

        %Save figure
        cd(fullfile(rns_results, 'Telemetry RNS Figures'))
        savefig(RNS_figure, [ID, ' Run ', num2str(run),' Channel ', num2str(channel),'.fig'])
    
    end
    
    %Plot additional points for 535BP
    if isequal(ID, '535BP') || isequal(ID, '536JD')
        
        scatter(new_telemetry_times,ones(1,size(new_telemetry_times,1))*-512, 'r*'); 
    
    end
    
    %Select telemetry times of interest (onset of final pulse)
    if isequal(ID, '471MH')
        
        %Remove first mark time in run for 471MH
        if ismember(run, [1,3,4,5,6,7,8,10,11,13,14,15])
           
           telemetry_times(1:2) = []; 
        
        %Run 9   
        elseif ismember(run, 9)
           
            telemetry_times(1:4) = []; 
            
        end
        
        %Store final telemetry times - Selecting the first deflection point
        %of the last of the four RNS markers 
        final_telemetry_times = telemetry_times(13:18:end);
               
    elseif isequal(ID, '535BP') || isequal(ID, '536JD')
        
        %Store final telemtry times
        final_telemetry_times = new_telemetry_times(19:24:end);
               
    elseif isequal(ID, '490KB')
        
        %Remove marks from run 1
        if ismember(run, 1)
            
            %Remove individual mark associated with missing mark set
            telemetry_times(409:414) = [];
            
        end
        
        %Remove marks from run 6
        if ismember(run, 6)
            
           %Add missing marks times
           telemetry_times = sort([telemetry_times', [37715,37716,91489,91490]])';
        
        end
        
        %Store final telemetry times
        final_telemetry_times = telemetry_times(19:24:end);
        
    else
        
        %Store final telemetry times
        final_telemetry_times = telemetry_times(19:24:end);
    
    end
    
    %% Find and store face times
    
    %Setup reject_idx
    reject_idx = zeros(length(run_pre_stim),1);
    
    %Remove first three trials with missing marks
    if ismember(run, 11) && isequal(ID, '471MH')
       
       %run_pre_stim(1:3) = [];
       %trial_ID(1:3) = [];
       reject_idx(1:3) = 1;
       NaN_telemetry_idx = nan(length(run_pre_stim),1);
       NaN_telemetry_idx(4:length(run_pre_stim),1) = final_telemetry_times;
       final_telemetry_times = NaN_telemetry_idx;
       
    %Remove first trial from first run because mark was not registered in RNS
    elseif ismember(run, 1) && isequal(ID, '535BP')
        
       %run_pre_stim(1) = [];
       %trial_ID(1) = [];
       reject_idx(1) = 1;
       NaN_telemetry_idx = nan(length(run_pre_stim),1);
       NaN_telemetry_idx(2:length(run_pre_stim),1) = final_telemetry_times;
       final_telemetry_times = NaN_telemetry_idx;

    %Remove trial corresponding to dropped mark trial   
    elseif ismember(run, 1) && isequal(ID, '490KB')
        
        %run_pre_stim(18) = [];
        %trial_ID(18) = [];
        reject_idx(18) = 1;
        NaN_telemetry_idx = nan(length(run_pre_stim),1);
        NaN_telemetry_idx(1:17,1) = final_telemetry_times(1:17,1);
        NaN_telemetry_idx(19:length(run_pre_stim),1) = final_telemetry_times(18:length(run_pre_stim)-1,1);
        final_telemetry_times = NaN_telemetry_idx;
        
    end
    
    %Check the number of telemetry times and pre-stim times are equal
    if not(isequal(length(final_telemetry_times), length(run_pre_stim)))
        
        error('Number of pre-stim times and telemetry times mismatch!')
        
    end
    
    %Add delays and prestim times to telemetry sample (not including
    %run_pre_stim that are member of the reject_idx)
    %face_times = final_telemetry_times + RNS_delay*sampling_rate + run_pre_stim(not(reject_idx))*sampling_rate;
    face_times = final_telemetry_times + RNS_delay*sampling_rate + run_pre_stim*sampling_rate;
        
    %% Cut CP and CnP epochs    
        
    %Initialize variable
    face_epochs = [];
    
    %Cut epochs
    for face = 1:length(face_times)    
        
        %loop over channels
        for channel = 1:4
            
            %Update reject idx 
            if isnan(final_telemetry_times(face))
                
                %Output variable - reject
                keep_or_reject = 'reject';
                
            else
                
                %Output variable - default is keep
                keep_or_reject = 'keep';
            
            end
            
            %If Subject ID 536JD run function below that eliminates bad
            %trials with epileptiform activity 
            if isequal(ID, '536JD')
               
               %Run keep or reject trial function
               [keep_or_reject] = remove_seizure_trials_function_536JD(run,face,channel);
                
            end
            
           %Skip rejected trials
           if isequal(keep_or_reject, 'reject')

                %continue
                
                %Make a NaN matrix to replace rejected trial values
                epoch = nan(1,1001);
           
           %Keep trial and cut epoch
           elseif isequal(keep_or_reject, 'keep')
               
                %Cut epoch for signal face event
                epoch = ECoG_data_combined{1,channel}(1,face_times(face)-(epoch_duration*sampling_rate):face_times(face)+(epoch_duration*sampling_rate));
                
           end
           
           %Aggregate epochs
           face_epochs(channel,:,face) = epoch;
           
        end
        
    end
    
    %Check the variables are equal sized
    if not(isequal(length(trial_ID),size(face_epochs,3),length(final_telemetry_times),length(face_times)))
        
        error('Variable sizes do not match!')
        
    end
    
    %Store the CP epochs
    all_CP_epochs = cat(3,all_CP_epochs, face_epochs(:,:,strcmp(trial_ID, 'Confirmed perceived')));
    all_CnP_epochs = cat(3,all_CnP_epochs, face_epochs(:,:,strcmp(trial_ID, 'Confirmed not perceived')));
    
end

%Remove bad CP and CnP trial from 490KB
if isequal(ID, '490KB')
    
   %Remove CP trials with massive artifact 
   %all_CP_epochs(:,:,30) = []; 
   all_CP_epochs(:,:,30) = NaN; 

   %Remove CnP trials with massive artifact 
   %all_CnP_epochs(:,:,[62,65,80,81]) = []; 
   all_CnP_epochs(:,:,[62,65,80,81]) = NaN; 
   
%Remove bad CP and CnP trial from 535BP   
elseif isequal(ID, '535BP')
    
    %Remove CP trials with massive artifact
    %all_CP_epochs(:,:,[54]) = [];
    all_CP_epochs(:,:,[54]) = NaN;
    
    %Remove CnP trials with massive artifact
    %all_CnP_epochs(:,:,[45]) = [];
    all_CnP_epochs(:,:,[45]) = NaN;
    
%Remove bad CP and CnP trial from 536JD  
elseif isequal(ID, '536JD')
    
    %Remove CnP trials with massive artifact
    %all_CnP_epochs(:,:,[67]) = [];
    all_CnP_epochs(:,:,[67]) = NaN;

end

%Save RNS epochs
cd(epoch_results)
save('RNS_CP_CnP_epochs.mat', 'all_CP_epochs', 'all_CnP_epochs', '-v7.3')

%{
%% Plot Mean CP and CnP Responses with EEG

%Load EEG Data
load(['P:\RNS Study\Analysis\Subject Analysis\',ID,'\Scalp EEG Analysis\CP_CnP_filtered_6000ms_epochs.mat']);

%Scalp channel to plot
scalp_channel = 19; %51

%Average over trials
mean_scalp_CP_epochs = mean(CP_epochs,3);
mean_scalp_CnP_epochs = mean(CnP_epochs,3);

%Resample from scalp data from 256 to 250
mean_channel_scalp_CP_epochs = resample(mean_scalp_CP_epochs(scalp_channel,:), 250,256);
mean_channel_scalp_CnP_epochs = resample(mean_scalp_CnP_epochs(scalp_channel,:), 250,256);

%Define time for plot
time = [-2000:1000/sampling_rate:2000];

%Loop over channels
for channel = 1:4
    
    %Calculate SEM
    CP_SEM = transpose(std(squeeze(all_CP_epochs(channel,:,:)),0, 2)/sqrt(size(all_CP_epochs,3)));
    CnP_SEM = transpose(std(squeeze(all_CnP_epochs(channel,:,:)),0, 2)/sqrt(size(all_CnP_epochs,3)));
    
    %Plot Mean RNS Epochs
    figure
    hold on
    ylim([-45 45])

    title([ID, ' RNS CP vs CnP - Channel #', num2str(channel)])
    xlabel('Time (ms)')
    ylabel('Voltage (uV)')
    
    %Plot traces
    RNS_CP = plot(time, nanmean(all_CP_epochs(channel,:,:),3), 'b');
    RNS_SEM_CP_1 = plot(time, nanmean(all_CP_epochs(channel,:,:),3)+ CP_SEM, '--b');
    RNS_SEM_CP_2 = plot(time, nanmean(all_CP_epochs(channel,:,:),3)- CP_SEM, '--b');
    
    RNS_CnP = plot(time, nanmean(all_CnP_epochs(channel,:,:),3), 'r');
    RNS_SEM_CnP_1 = plot(time, nanmean(all_CnP_epochs(channel,:,:),3)+ CnP_SEM, '--r');
    RNS_SEM_CnP_2 = plot(time, nanmean(all_CnP_epochs(channel,:,:),3)- CnP_SEM, '--r');
    
    EEG_CP = plot(time, mean_channel_scalp_CP_epochs(250:1250), 'g');
    EEG_CnP = plot(time, mean_channel_scalp_CnP_epochs(250:1250), 'y');
    
    %Reference line
    plot([0,0],[-45,45],'k')
    
    %Generate legend
    legend([RNS_CP, [],[],RNS_CnP, [],[], EEG_CP, EEG_CnP,[]], 'RNS CP', 'RNS CnP', 'EEG CP', 'EEG CnP')
    
    %Save figure
    cd(rns_results)
    savefig(['RNS_EEG_ERPs_CP_CnP_channel_',num2str(channel)])
    
end
%}