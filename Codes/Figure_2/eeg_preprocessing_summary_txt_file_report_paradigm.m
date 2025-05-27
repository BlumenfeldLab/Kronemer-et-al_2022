%% EEG Preprocessing Summary Output

%This code will create a text files that summarizes the preprocessing
%rejection stats - session level and epoch level rejections

%Written by: Sharif I. Kronemer
%Date: 1/28/2021
%Modified: 4/29/2021

function eeg_preprocessing_summary_txt_file_report_paradigm(ID, EEG_data, bad_channel_dir, preprocessed_dir)

%Load bad channel/sample/epoch info
cd(bad_channel_dir)
load('Session_bad_channels_samples.mat')
load('Main_events_bad_epochs_samples.mat')

%Create text file 
cd(preprocessed_dir)
test_file = fopen('EEG_preprocessing_summary_stats.txt','wt');
fprintf(test_file,'%s\n\r\n',['EEG Preprocessing Summary: ', ID]);

%% Session-level Rejection Summary

fprintf(test_file,'%s\r\n',['**Session Rejections**']);

%Loop through the raw_files that equals the number of available behavioral run csv files
for raw_files = 1:size(EEG_data,2)
    
    %Rename variables
    eval(['bad_channels = Session_',num2str(raw_files),'_bad_channels;'])
    eval(['bad_samples = Session_',num2str(raw_files),'_bad_samples;'])

    %Print session rejection results
    fprintf(test_file,'\n%s\r\n',['Session ', num2str(raw_files),' - Channels rejected (%)']);
    fprintf(test_file,'%f\n',(length(bad_channels)/256)*100);
    fprintf(test_file,'\r\n%s\r\n',['Session ', num2str(raw_files), ' - Samples rejected (%)']);
    fprintf(test_file,'%f\n', (1-(sum(bad_samples)/length(bad_samples)))*100);

end

%% Epoch-level Rejection Summary

fprintf(test_file,'\n%s\r\n',['**Epoch Rejections**']);

%Print epoch rejection results
fprintf(test_file,'\n%s\r\n',['Face epochs rejected (%)']);
fprintf(test_file,'%f\n',(sum(All_faces_bad_epochs_idx)/length(All_faces_bad_epochs_idx))*100);
% fprintf(test_file,'\r\n%s\r\n',['Button press epochs rejected (%)']);
% fprintf(test_file,'%f\n',(sum(All_button_presses_bad_epochs_idx)/length(All_button_presses_bad_epochs_idx))*100);
% fprintf(test_file,'\r\n%s\r\n',['Question epochs rejected (%)']);
% fprintf(test_file,'%f\n',(sum(All_questions_bad_epochs_idx)/length(All_questions_bad_epochs_idx))*100);

end