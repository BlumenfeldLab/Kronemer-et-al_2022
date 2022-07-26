%% No Report Paradigm Relevant Dataset Classifier Training

%This script loads relevant datasets from no report dataset (CP/CnP trials)
%and trainings a stacked ensemble classifier that can be used to test on
%irrelevant dataset

%General architecture
%(1) Train base models (Gradient boost, Linear SVM, and Gaussian/RBF SVM)
% individually for each data type (interpolated pupil diameter, tsfresh
% features from interpolated pupil diameter after applying t-test feature
% selection, raw Y-axis gaze, raw X-axis gaze, binary (1/0) blink event
% logical, and binary (1/0) microsaccade event logical)
%(2)Stacked model is a linear SVM that is trained on the confidence scores
%emerging from each base model for each data type (3 base models x 6 data
%types) 

%Written by: Sharif I. Kronemer
%Date: 4/18/2021
%Modified: 7/5/2022

clear

%% Prompts

%Select run location
prompt_1 = 'Running code local or server [l, s]: ';
run_location = input(prompt_1,'s');

%Relevant condition
prompt_2 = 'Center and/or quadrant relevant [c, q, combine]: ';
condition = input(prompt_2,'s');

%Center relevant
if isequal(condition, 'c')

    relevant_name = 'Center Relevant';
    relevant_location = {'Center Relevant'};

%Quadrant relevant
elseif isequal(condition, 'q')
    
    relevant_name = 'Quadrant Relevant';
    relevant_location = {'Quadrant Relevant'};
    
%Combine Center and Quadrant Relevant   
elseif isequal(condition, 'combine')
    
    relevant_name = 'Center and Quadrant Relevant';
    relevant_location = {'Center Relevant', 'Quadrant Relevant'};
    
end

%Modality
prompt_3 = 'Modality [EEG, MRI]: ';
modality = input(prompt_3,'s');

%Opacity
prompt_4 = 'Stimulus opacity [blank,threshold]: ';
opacity = input(prompt_4,'s');

%%  Define directories and paths

%Directory save filename 
filename = 'Base_boosted_tree_lin_gaus_SVM_stacked_lin_SVM';
%filename = 'Base_boosted_tree_lin_gaus_SVM_stacked_lin_SVM_no_tsfresh'; % Train/Test model without tsfresh features

%Save direcotry
if isequal(run_location,'l')
    

elseif isequal(run_location,'s')
    
    %Add paths
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/Behavioral Analysis')
    addpath(genpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/No Report Paradigm Machine Learning/Training Testing and Plotting Functions'))

    %Save directory
    save_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/Machine Learning Analysis',[modality,' Session'],relevant_name,filename);

    %Subject EEG and MRI dir
    EEG_sub_dir = '/mnt/Data28/HNCT No Report Paradigm/Subject Analysis EEG';
    MRI_sub_dir = '/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI';

end
    
%Make save directory if it does not exist
if not(exist(save_dir))
    
    mkdir(save_dir)
    
else
    
    warning('***Save directory already exists. Change filename if you do not intend to overwrite previous results.***')

end

%% STEP 1: Load EEG and Eyelink data - Combine over subjects

%Loop over relevant locations
for n = 1:length(relevant_location)
    
    %Define current location
    current_location = relevant_location{n};

    %Load subject EEG and EEG-eyelink data and add to a group variable. Output
    %the group extracted eyelink, EEG, and true CP/CnP labels. 
    if isequal(modality, 'EEG')

        disp(['Extract EEG Eyelink Data - ',current_location])

        %Extract the relevant eyelink trials - EEG
        [report_data_pupil, report_data_tsfreshpupil, report_data_Ygaze, report_data_Xgaze, report_data_blink, report_data_microsac, ...
            true_labels, stimulus_opacity, group_subject, group_EEG_epoch_count, number_of_EEG_subjects, subject_trial_index] = ...
            transfer_combine_eyelink_EEG_relevant_epochs_function(EEG_sub_dir, current_location, run_location);

    elseif isequal(modality, 'MRI')

        disp(['Extract MRI Eyelink Data - ',current_location])

        %Extract the relevant eyelink trials  - MRI
        [report_data_pupil, report_data_tsfreshpupil, report_data_Ygaze, report_data_Xgaze, report_data_blink, report_data_microsac, true_labels, ...
            group_subject, group_run_number, group_trial_number, group_condition, group_relevant_location, group_BOLD_epoch_count, number_of_MRI_subjects, subject_trial_index] = ...
            combine_eyelink_MRI_relevant_epochs_function(MRI_sub_dir, current_location, run_location, opacity);

    elseif isequal(modality, 'combine')

        %NOTE: Optimization testing of the classifier performance suggests the
        %current dataset is best predicted when EEG and MRI sessions are
        %treated separately

        disp('Extract and Combine EEG and MRI Eyelink Data')

        %Extract the relevant eyelink trials - EEG
        [report_data_pupil, report_data_tsfreshpupil, report_data_Ygaze, report_data_Xgaze, report_data_blink, report_data_microsac, ...
            true_labels, stimulus_opacity, group_subject, group_EEG_epoch_count, number_of_EEG_subjects] = combine_eyelink_EEG_relevant_epochs_function(EEG_sub_dir, current_location, run_location);

        %Rename EEG eyelink variables
        report_data_pupil_eeg = report_data_pupil;
        report_data_tsfreshpupil_eeg = report_data_tsfreshpupil;
        report_data_Ygaze_eeg = report_data_Ygaze;
        report_data_Xgaze_eeg = report_data_Xgaze;
        report_data_blink_eeg = report_data_blink;
        report_data_microsac_eeg = report_data_microsac;
        true_labels_eeg = true_labels;

        %Extract the relevant eyelink trials  - MRI
        [report_data_pupil, report_data_tsfreshpupil, report_data_Ygaze, report_data_Xgaze, report_data_blink, report_data_microsac, true_labels, ...
            group_subject, group_run_number, group_trial_number, group_condition, group_relevant_location, group_BOLD_epoch_count, number_of_MRI_subjects] = ...
            combine_eyelink_MRI_relevant_epochs_function(MRI_sub_dir, current_location, run_location);

        %Combine EEG variables to MRI variables
        report_data_pupil = cat(1, report_data_pupil, report_data_pupil_eeg);
        report_data_tsfreshpupil = cat(1, report_data_tsfreshpupil, report_data_tsfreshpupil_eeg);
        report_data_Ygaze = cat(1, report_data_Ygaze, report_data_Ygaze_eeg);
        report_data_Xgaze = cat(1, report_data_Xgaze, report_data_Xgaze_eeg);    
        report_data_blink = cat(1, report_data_blink, report_data_blink_eeg);
        report_data_microsac = cat(1, report_data_microsac, report_data_microsac_eeg);
        true_labels = [true_labels; true_labels_eeg];

    end  
    
    %Combine data across locations
    if n == 1 && length(relevant_location) == 2
        
        report_data_pupil_loc1 = report_data_pupil;
        report_data_tsfreshpupil_loc1 = report_data_tsfreshpupil;
        report_data_Ygaze_loc1 = report_data_Ygaze;
        report_data_Xgaze_loc1 = report_data_Xgaze;    
        report_data_blink_loc1 = report_data_blink;
        report_data_microsac_loc1 = report_data_microsac;
        true_labels_loc1 = true_labels;
       
    elseif n == 2 && length(relevant_location) == 2
        
        report_data_pupil = cat(1, report_data_pupil, report_data_pupil_loc1);
        report_data_tsfreshpupil = cat(1, report_data_tsfreshpupil, report_data_tsfreshpupil_loc1);
        report_data_Ygaze = cat(1, report_data_Ygaze, report_data_Ygaze_loc1);
        report_data_Xgaze = cat(1, report_data_Xgaze, report_data_Xgaze_loc1);    
        report_data_blink = cat(1, report_data_blink, report_data_blink_loc1);
        report_data_microsac = cat(1, report_data_microsac, report_data_microsac_loc1);
        true_labels = [true_labels; true_labels_loc1];
        
    end
    
    
end

%% Preprocess Training Data

%Define data epoch (3s post stimulus)
data_epoch = 6001:9000;

%Down sample - number of samples 1000ms to 50 samples
imresize_num = 150;

%Cut the epoch range
report_data_pupil_cut = report_data_pupil(:,data_epoch);
report_data_Ygaze_cut = report_data_Ygaze(:,data_epoch);
report_data_Xgaze_cut = report_data_Xgaze(:,data_epoch);
report_data_blink_cut = report_data_blink(:,data_epoch);
report_data_microsac_cut = report_data_microsac(:,data_epoch);

%Downsample the data for computational efficiency
report_data_pupil_cut = imresize(report_data_pupil_cut, [size(report_data_pupil_cut,1) imresize_num]); 
report_data_Ygaze_cut = imresize(report_data_Ygaze_cut, [size(report_data_Ygaze_cut,1) imresize_num]); 
report_data_Xgaze_cut = imresize(report_data_Xgaze_cut, [size(report_data_Xgaze_cut,1) imresize_num]); 
report_data_blink_cut = imresize(report_data_blink_cut, [size(report_data_blink_cut,1) imresize_num]);  
report_data_microsac_cut = imresize(report_data_microsac_cut, [size(report_data_microsac_cut,1) imresize_num]);  

%Return blink and microsac to binary after imresize
report_data_blink_cut(report_data_blink_cut ~= -1) = 1;
report_data_microsac_cut(report_data_microsac_cut ~= 0) = 1;

%% Feature Selection

disp('Running feature selection')

%Feature selection alpha value
alpha_value = 0.05;

%Pupil tsfresh feature selection
[h,p,ci,stat] = ttest2(report_data_tsfreshpupil(find(true_labels==1),:),report_data_tsfreshpupil(find(true_labels==2),:),'alpha',alpha_value);    
report_data_tsfreshpupil_cut = report_data_tsfreshpupil(:,find(h==1));

%Store significant pupilfeatures
selected_tsfreshpupil_features = find(h==1);

%% Train Base/Level 0 Models

%Create data name variable 
%data_names = {'pupil','Ygaze','Xgaze','blink','microsac'};  %No tsfresh
data_names = {'pupil','tsfreshpupil','Ygaze','Xgaze','blink','microsac'}; %With tsfresh

%Model type cell
model_names = {'boosted_tree','SVM_gaus','SVM_lin'};

%Initialize variable
base_model_scores = [];
base_model_labels = [];
base_model_names = [];

%Setup data partitioning to be used for cross validation
rng('default') % For reproducibility 
cv = cvpartition(true_labels,"KFold",10); 

%Loop over data
for types = 1:length(data_names)

    tic
    
    %Rename data variable
    current_data_name = data_names{types};
    
    disp(['Running ',current_data_name])
    
    %Initialize training data
    training_data = [];

    %Rename current variable to generic name for training
    eval(['training_data = report_data_',current_data_name,'_cut;']);
  
    %Train hetergeneous models - base models
    disp('Train base models')
   
    %Boosted ensemble deceision tree
    classModels{1} = fitcensemble(training_data, true_labels);
    
    %Gaussian SVM
    classModels{2} = fitcsvm(training_data, true_labels, 'Standardize', true, 'KernelFunction', 'gaussian');
    
    %Linear SVM 
    classModels{3} = fitcsvm(training_data, true_labels, 'Standardize', true, 'KernelFunction', 'linear');
    
    %Loop over models and store scores
    for model = 1:numel(classModels)

        disp(['Cross validate model ',num2str(model)])
        
        %Select current model among the base models
        current_model = model_names{model};
        
        %Save data/model type name 
        base_model_names{end+1,1} = [current_data_name,' ',current_model];
        
        %Cross-validate
        current_CV_Model = crossval(classModels{model},'cvpartition',cv); 
        
        %Predict classes and generate base model confidence scores
        [current_labels, current_scores] = kfoldPredict(current_CV_Model);

        %Store base model confidence scores and labels
        base_model_scores = [base_model_scores, current_scores(:,1)]; 
        base_model_labels = [base_model_labels, current_labels]; 
        
        %Rename cross-validated base model
        eval(['CV_model_',current_model,'_',current_data_name, ' = current_CV_Model;'])
        
    end
    
    toc

end

% Train Stacked/Level 1 Model

%Linear SVM 
stacked_model = fitcsvm(base_model_scores, true_labels, 'Standardize', true);

%Cross-validate
CV_stacked_model = crossval(stacked_model,'cvpartition',cv);

%Predict classes and generate stacked model confidence scores
[stacked_labels, stacked_scores] = kfoldPredict(CV_stacked_model);

%Save model variables
cd(save_dir)
save('Classifier_models_scores_labels.mat', 'CV_model*', 'classModels','CV_stacked_model','base_model_scores', 'base_model_labels' ,'base_model_names', 'stacked_labels' ,'stacked_scores', 'model_names', ...
    'data_names', 'true_labels' ,'group_subject', 'alpha_value', 'imresize_num' ,'data_epoch' ,'selected_tsfreshpupil_features' ,'modality' ,'relevant_location','-v7.3');

% Plot Stacked Model Confidence Score Histogram

%Setup figure
figure
hold on

title(['Stacked SVM Model - ', modality, ' - ', relevant_name], 'Interpreter', 'none')
xlabel('Stacked SVM Scores')
ylabel('Number of Trials')
xlim([-8,8])

%Plot histogram of true CP and CnP trial confidence scores
histogram(stacked_scores(true_labels == 1),50,'FaceColor','b')
histogram(stacked_scores(true_labels == 2),50,'FaceColor','r')

%Save figure
cd(save_dir)
savefig(['Stacked_SVM_model_',modality,'_confidence_score_histogram.fig'])

close

% Plot Stacked Model PPV/NPV vs Retention

%Confidence Score Thresholds
confidence_score_thresholds = [0:0.25:3];    

%PPV and Proportion Predicted 
PPV = figure
hold on

title(['Confirmed Perceived and Not Perceived - ', relevant_name,' - ',modality], 'Interpreter', 'none')

xlabel('PPV / NPV')
ylabel('Proportion of Trials')
xlim([0.65 1])
ylim([0 1.2])

%Plot reference lines
refline_1 = plot([0.55 1],[0.5 0.5],'Color','k','LineStyle','--')
refline_2 = plot([0.9 0.9],[0 1.5],'Color','k','LineStyle','--')

%Loop over confidence scores
for thres = 1:length(confidence_score_thresholds)

    %Define threshold 
    current_threshold = confidence_score_thresholds(thres);
    
    %Confirmed Perceived

    %Correctly predicted class / all the predicted class
    PPV_stacked = length(find(stacked_labels(stacked_scores(:,1) > current_threshold) == true_labels(stacked_scores(:,1) > current_threshold)...
            & stacked_labels(stacked_scores(:,1) > current_threshold) == 1))/length(find(stacked_labels(stacked_scores(:,1) > current_threshold) == 1));

    %Find number of predicted trials
    number_predicted_trials_stacked = length(find(stacked_labels(stacked_scores(:,1) > current_threshold) == 1)); 

    %Percentage of predicted trials
    CP_retention_stacked = number_predicted_trials_stacked/length(find(true_labels == 1));  

    %PPV vs retention
    class_1 = scatter(PPV_stacked,CP_retention_stacked,'filled', 'k')
    
    %Confirmed not Perceived

    %Correctly predicted class / all the predicted class
    NPV_stacked = length(find(stacked_labels(stacked_scores(:,1) < -current_threshold) == true_labels(stacked_scores(:,1) < -current_threshold)...
            & stacked_labels(stacked_scores(:,1) < -current_threshold) == 2))/length(find(stacked_labels(stacked_scores(:,1) < -current_threshold) == 2));

    %Find number of predicted trials
    number_predicted_trials_stacked = length(find(stacked_labels(stacked_scores(:,1) < -current_threshold) == 2)); 

    %Percentage of predicted trials
    CnP_retention_stacked = number_predicted_trials_stacked/length(find(true_labels == 2));  

    %NPV vs retention
    class_2 = scatter(NPV_stacked,CnP_retention_stacked,'filled', 'r')

    %Setup legend
    legend_cell{thres} =['Score thres ',num2str(current_threshold)];
    
    % accuracy
    disp(sum(true_labels==stacked_labels)/length(true_labels))

end

%Plot legend and figure text
legend([class_1,class_2],{'Confirmed perceived (PPV)','Confirmed not perceived (NPV)'})
text(0.66,0.2, {'Stacked SVM Score Thresholds:','0.0   0.25   0.5','0.75   1.0   1.25',...
    '1.5   1.75   2.0','2.25   2.5   2.75   3.0'})
 
%Save figure
cd(save_dir)
savefig(['PPV_NPV_vs_retention_stacked_SVM_',modality,'_figure.fig'])


%% Plot ROC 1-Specificity vs Sensitivity

%ROC Plot
ROC = figure
hold on

title(['ROC Base and Stacked Models - ', relevant_name,' - ',modality], 'Interpreter', 'none')

xlabel('1-Specificity')
ylabel('Sensitivity')
xlim([0 1])
ylim([0 1])

%Plot 45deg line
plot([0 1],[0 1],'k')

%Loop over all base modeles
for model = 1:size(base_model_scores,2)
    
    %Find the ROC values; Calculate ROC by Matlab Function that uses cross validation
    [x,y] = perfcurve(true_labels,base_model_scores(:,model),1);
    
    %Decide on the line color
    
    %Pupil 
    if contains(base_model_names{model},'pupil ') & not(contains(base_model_names{model},'tsfresh'))
        
        %Plot current base model ROC values
        pupilbase_plot = plot(x,y,'c')
    
    %TsFresh Pupil
    elseif contains(base_model_names{model},'tsfresh')
        
        %Plot current base model ROC values
        tsfreshbase_plot = plot(x,y,'g')
        
    %Gaze data X and Y
    elseif contains(base_model_names{model},'gaze')
        
        %Plot current base model ROC values
        gazebase_plot = plot(x,y,'m')
        
    %Blink data
    elseif contains(base_model_names{model},'blink')
        
        %Plot current base model ROC values
        blinkbase_plot = plot(x,y,'b')
        
    %Microsaccade data
    elseif contains(base_model_names{model},'microsac')
        
        %Plot current base model ROC values
        microsacbase_plot = plot(x,y,'y')
        
    end

end

%Find the ROC values for the stacked/ensemble score
[x,y,T,AUC,OPTROCPT] = perfcurve(true_labels,stacked_scores(:,1),1);

%Plot stacked ROC values
stacked_plot = plot(x,y,'r')

%Plot the optimal point of the ROC 
opt_plot = plot(OPTROCPT(1),OPTROCPT(2),'ko')

legend([pupilbase_plot, tsfreshbase_plot, gazebase_plot, blinkbase_plot, stacked_plot],...
    {'Pupil Base Model', 'Tsfresh Base Model', 'Gaze Base Model', 'Blink Base Model',...
    ['Stacked Ensemble Model - AUC = ', num2str(AUC)]},'Location','best')

%{
%Uses two confidence threshold boundaries

%Confidence Score Thresholds
confidence_score_thresholds = [0:0.25:2];    

%Loop over confidence scores
for thres = 1:length(confidence_score_thresholds)

    %Define threshold 
    current_threshold = confidence_score_thresholds(thres);    
    
    %TP
    TP = length(find(stacked_labels(stacked_scores(:,1) > current_threshold) == true_labels(stacked_scores(:,1) > current_threshold)...
            & stacked_labels(stacked_scores(:,1) > current_threshold) == 1));
     
    %TN
    TN = length(find(stacked_labels(stacked_scores(:,1) < -current_threshold) == true_labels(stacked_scores(:,1) < -current_threshold)...
            & stacked_labels(stacked_scores(:,1) < -current_threshold) == 2));    
    
    %FP
    FP = length(find(stacked_labels(stacked_scores(:,1) > current_threshold) == 1 & true_labels(stacked_scores(:,1) > current_threshold) == 2));
    
    %FN
    FN = length(find(stacked_labels(stacked_scores(:,1) < -current_threshold) == 2 & true_labels(stacked_scores(:,1) < -current_threshold) == 1));
        
    %Senstivity = TP/TP+FN (NOTE: PPV is different = TP/TP+FP)
    sensitivity_stacked = TP/(TP+FN);
    
    %Specificity = TN/TN+FP
    specificity_stacked = TN/(TN+FP);  
    
    %1-Specificity vs Sensitivity
    ROC_point = scatter(1-specificity_stacked,sensitivity_stacked,'filled', 'k')

    %Setup legend
    legend_cell{thres} =['Score thres ',num2str(current_threshold)];
    
end
%}

%{
%Uses one confidence threshold boundary

confidence_score_thresholds = [-3:0.25:3];    

%Loop over confidence scores
for thres = 1:length(confidence_score_thresholds)

    %Define threshold 
    current_threshold = confidence_score_thresholds(thres);    
    
    %TP
    TP = length(find(stacked_labels(stacked_scores(:,1) > current_threshold) == true_labels(stacked_scores(:,1) > current_threshold)...
            & stacked_labels(stacked_scores(:,1) > current_threshold) == 1));
     
    %TN
    TN = length(find(stacked_labels(stacked_scores(:,1) < current_threshold) == true_labels(stacked_scores(:,1) < current_threshold)...
            & stacked_labels(stacked_scores(:,1) < current_threshold) == 2));    
    
    %FP
    FP = length(find(stacked_labels(stacked_scores(:,1) > current_threshold) == 1 & true_labels(stacked_scores(:,1) > current_threshold) == 2));
    
    %FN
    FN = length(find(stacked_labels(stacked_scores(:,1) < current_threshold) == 2 & true_labels(stacked_scores(:,1) < current_threshold) == 1));
        
    %Senstivity = TP/TP+FN (NOTE: PPV is different = TP/TP+FP)
    sensitivity_stacked = TP/(TP+FN);
    
    %Specificity = TN/TN+FP
    specificity_stacked = TN/(TN+FP);  
    
    %1-Specificity vs Sensitivity
    ROC_point = scatter(1-specificity_stacked,sensitivity_stacked,'filled', 'k')

    %Setup legend
    legend_cell{thres} =['Score thres ',num2str(current_threshold)];
    
end
%}

%Save figure
cd(save_dir)
savefig(['ROC_base_stacked_models_',modality,'_figure.fig'])

close

%% Plot Eyelink Timecourses - Predicted vs True Class

%Exclude trials with blink at stimulus onset? Yes = 1; No = 0
exclude_blink_trials = 1;

%Run plotting function
eyelink_timecourse_predicted_relevant_stimuli_function(save_dir, relevant_name, modality, report_data_pupil, report_data_blink, report_data_microsac, exclude_blink_trials)
