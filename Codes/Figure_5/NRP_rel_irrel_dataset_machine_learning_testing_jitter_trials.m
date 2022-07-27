%% No Report Paradigm Relevant Irrelevant Dataset Classifier Testing - Post-Stimulus Jitter Trials

% Jitter epochs predicted perceived/not perceived classification code
% variant 

%NOTE: Due to processing pipeline differences for EEG Irrelevant and
%Relevant Center and Quadrant conditions are considered separately. But,
%for MRI only relevant conditions are considered in which the irrelevant
%class for that relevant condition is automatically included. 

%Written by: Sharif I. Kronemer
%Date: 4/18/2021
%Modified: 6/5/2022

clear

%% Prompts

%Select run location
prompt_1 = 'Running code local or server [l, s]: ';
run_location = input(prompt_1,'s');

%Modality
prompt_3 = 'Modality [EEG, MRI]: ';
modality = input(prompt_3,'s');

%If EEG condition
if isequal(modality, 'EEG')
    
    rel_condition = input('Relevant or Irrelevant [r, i]: ','s');
    
else
    
    rel_condition = 'r';

end

%Location condition
prompt_2 = 'Center or Quadrant [c, q]: ';
condition = input(prompt_2,'s');

%Center relevant
if isequal(condition, 'c')

    irrelevant_location = 'Center Irrelevant';
    relevant_location = 'Center Relevant';

%Quadrant relevant
elseif isequal(condition, 'q')
    
    irrelevant_location = 'Quadrant Irrelevant';
    relevant_location = 'Quadrant Relevant';
    
end



%%  Define directories and paths

%Directory save filename 
model_filename = 'Base_boosted_tree_lin_gaus_SVM_stacked_lin_SVM_no_tsfresh'; %'Base_boosted_tree_lin_gaus_SVM_stacked_lin_SVM';
save_filename = 'Base_boosted_tree_lin_gaus_SVM_stacked_lin_SVM_jitter_new'; 

%Report data filename
report_data_filename = 'Group_eyelink_CP_CnP_cent_quad_threshold_EEG_MRI_data.mat';

%Save direcotry
if isequal(run_location,'l')

elseif isequal(run_location,'s')
    
    %Add paths
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/Machine Learning/Ensemble SVM Approach')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/Machine Learning')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/Behavioral Analysis')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/No Report Paradigm Machine Learning/Training Testing and Plotting Functions/MRI functions')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/No Report Paradigm Machine Learning/Training Testing and Plotting Functions/EEG functions')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/No Report Paradigm Machine Learning/Training Testing and Plotting Functions/Eyelink functions')
    
    %Trained classifier directory
    %model_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/Machine Learning Analysis',[modality,' Session'],relevant_location,model_filename);
    model_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/Machine Learning Analysis',[modality,' Session'],'Center and Quadrant Relevant',model_filename);

    if strcmp(modality, 'EEG')
        
        %Save directory
        if strcmp(rel_condition, 'i')

            save_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/Machine Learning Analysis',[modality,' Session'],irrelevant_location,save_filename);

        elseif strcmp(rel_condition, 'r')

            save_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/Machine Learning Analysis',[modality,' Session'],relevant_location,save_filename);

        end
    
    else
        
        save_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/Machine Learning Analysis',[modality,' Session'],relevant_location,save_filename);
       
    end
    
    %Report group Eyelink data
    report_data_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Data/EyeLink';
    
    %Subject dir
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

%Load subject EEG and EEG-eyelink data and add to a group variable. Output
%the group extracted eyelink, EEG, and true CP/CnP labels. 
if isequal(modality, 'EEG')
    
    disp('Extract EEG Eyelink Data')
    
    % Task Irrelevant
    if strcmp(rel_condition,'i')

        %Extract the irrelevant eyelink trials - EEG
        [noreport_data_pupil, noreport_data_tsfreshpupil, noreport_data_Ygaze, noreport_data_Xgaze, noreport_data_blink, noreport_data_microsac,...
            group_subject, group_run_number, group_trial_number, group_condition, group_relevant_location, group_relevant_trial, number_of_EEG_subjects, subject_trial_opacity, subject_trial_index] = ...
            combine_eyelink_EEG_irrelevant_epochs_jitter_function(EEG_sub_dir, irrelevant_location, run_location);

    % Task Relevant
    elseif strcmp(rel_condition,'r')

        [noreport_data_pupil, noreport_data_tsfreshpupil, noreport_data_Ygaze, noreport_data_Xgaze, noreport_data_blink, noreport_data_microsac, true_labels, ...
            group_subject, group_run_number, group_trial_number, group_condition, group_relevant_location, group_relevant_trial, number_of_EEG_subjects, subject_trial_opacity, subject_trial_index] = ...
            combine_eyelink_EEG_relevant_epochs_jitter_function(EEG_sub_dir, relevant_location, run_location);

    end
   
elseif isequal(modality, 'MRI')
    
    disp('Extract MRI Eyelink Data')

    %Extract the irrelevant eyelink trials - MRI

    [noreport_data_pupil, noreport_data_Ygaze, noreport_data_Xgaze, noreport_data_blink, noreport_data_microsac, group_subject, ...
        group_run_number, group_trial_number, group_condition, group_relevant_location, group_relevant_trial, group_BOLD_epoch_count, ...
        number_of_MRI_subjects, subject_trial_index, subject_trial_opacity] = combine_eyelink_MRI_rel_irrel_epochs_jitter_function(MRI_sub_dir, relevant_location, run_location);

end   

%% Preprocess Training Data

%Load trained model
cd(model_dir)
load Classifier_models_scores_labels.mat CV* selected_tsfreshpupil_features

%Define data epoch (3s post stimulus)
data_epoch = 9001:12000;

%Down sample - number of samples 1000ms to 50 samples
imresize_num = 150;

%Cut the epoch range
noreport_data_pupil_cut = noreport_data_pupil(:,data_epoch);
noreport_data_Ygaze_cut = noreport_data_Ygaze(:,data_epoch);
noreport_data_Xgaze_cut = noreport_data_Xgaze(:,data_epoch);
noreport_data_blink_cut = noreport_data_blink(:,data_epoch);
noreport_data_microsac_cut = noreport_data_microsac(:,data_epoch);

%Downsample the data so the code runs faster (every 180 samples)
noreport_data_pupil_cut = imresize(noreport_data_pupil_cut, [size(noreport_data_pupil_cut,1) imresize_num]); 
noreport_data_Ygaze_cut = imresize(noreport_data_Ygaze_cut, [size(noreport_data_Ygaze_cut,1) imresize_num]); 
noreport_data_Xgaze_cut = imresize(noreport_data_Xgaze_cut, [size(noreport_data_Xgaze_cut,1) imresize_num]); 
noreport_data_blink_cut = imresize(noreport_data_blink_cut, [size(noreport_data_blink_cut,1) imresize_num]);  
noreport_data_microsac_cut = imresize(noreport_data_microsac_cut, [size(noreport_data_microsac_cut,1) imresize_num]);  

%Return blink and microsac to binary after imresize
noreport_data_blink_cut(noreport_data_blink_cut ~= 0) = 1;
noreport_data_microsac_cut(noreport_data_microsac_cut ~= 0) = 1;

% Feature Selection according to relevant/labeled data tsfresh pupil feature selection
%noreport_data_tsfreshpupil_cut = noreport_data_tsfreshpupil(:,selected_tsfreshpupil_features);

%Replace tsfresh features with zeros
%noreport_data_tsfreshpupil_cut = zeros(size(noreport_data_tsfreshpupil_cut,1),size(noreport_data_tsfreshpupil_cut,2));

%% Train Base/Level 0 Models

%Data types (Note: order of data types matter to match with the models)
data_names =  {'pupil','Ygaze','Xgaze','blink','microsac'}; 

%Model types
model_names = {'boosted_tree','SVM_gaus','SVM_lin'};

%Initialize variable
base_model_scores = [];
base_model_names = [];

%Loop over data types
for types = 1:length(data_names)
    
    tic

    %Rename data variable
    current_data_name = data_names{types};
    
    disp(['Running ',current_data_name])
    
    %Initialize training data
    training_data = [];

    %Rename current variable to generic name for training
    eval(['training_data = noreport_data_',current_data_name,'_cut;']);
    
    %Loop over model types
    for model = 1:length(model_names)
        
        %Save data/model type name 
        base_model_names{end+1,1} = [current_data_name,' ',model_names{model}];

        %Rename particular trained model/data name
        eval(['CV_model = CV_model_',model_names{model},'_',current_data_name,';'])

        %Initialize variable
        all_models_labels = [];
        all_models_scores = [];
        consensus_labels = [];

        %Loop through all the models generated with cross validation above
        for model_num = 1:10

            %Select a fold within the ensemble model
            compacted_CV_Model = CV_model.Trained{model_num};

            %Test model on no report data - output label and score variables
            [label,score] = predict(compacted_CV_Model, training_data);

            %Store all the lables and scores for each fold
            all_models_labels = [all_models_labels, label];
            all_models_scores = [all_models_scores, score(:,1)];

        end 

        %Find prediction consensus across model predictions 
        %Note: Uses majority rule over models to select predicted class
        %label; if there is equal predictions between class 1 and 2
        %predicted class is 2
        for trial = 1:size(all_models_labels,1)

            %If the majority of models select class 1
            if length(find(all_models_labels(trial,:) == 1)) > 5

                %Consensus is class 1
                consensus_labels(trial,1) = 1;

            else

                %Consensus is class 2
                consensus_labels(trial,1) = 2;

            end

        end
        
        %Rename consensus labels
        eval(['consensus_labels_',model_names{model},'_',current_data_name,' = consensus_labels;'])
        
        %Average scores across models
        avg_models_score = mean(all_models_scores,2);

        %Aggregate mean scores over all models/data
        base_model_scores = [base_model_scores, avg_models_score];

    end
    
    toc

end

%% Train Stacked/Level 1 Model

%Initialize variable
stacked_SVM_all_models_labels = [];
stacked_SVM_all_models_scores = [];
stacked_SVM_consensus_labels = [];

%Loop through all the models generated with cross validation above
for model_num = 1:10

    %Select a fold within the ensemble model
    compacted_CV_Model = CV_stacked_model.Trained{model_num};

    %Test model on no report data - output label and score variables
    [label,score] = predict(compacted_CV_Model, base_model_scores);

    %Store all the lables and scores for each fold
    stacked_SVM_all_models_labels = [stacked_SVM_all_models_labels, label];
    stacked_SVM_all_models_scores = [stacked_SVM_all_models_scores, score(:,1)];

end 

%Find prediction consensus across model predictions 
%Note: Uses majority rule over models to select predicted class label;
%if there is equal predictions between class 1 and 2 predicted class is 2
for trial = 1:size(stacked_SVM_all_models_labels,1)

    %If the majority of models select class 1
    if length(find(stacked_SVM_all_models_labels(trial,:) == 1)) > 5

        %Consensus is class 1
        stacked_SVM_consensus_labels(trial,1) = 1;

    else

        %Consensus is class 2
        stacked_SVM_consensus_labels(trial,1) = 2;

    end

end

%Average scores across models
stacked_SVM_avg_scores = mean(stacked_SVM_all_models_scores,2);

%Save model variables
cd(save_dir)

if isequal(modality,'EEG')
           
    % Task Irrelevant
    if strcmp(rel_condition,'i')
        
        save Classifier_scores_labels_jitter.mat subject_trial_index stacked_SVM_avg_scores stacked_SVM_all_models_scores ...
            stacked_SVM_all_models_labels stacked_SVM_consensus_labels base_model_scores base_model_names ...
            imresize_num data_epoch modality irrelevant_location consensus_labels_* group* subject_trial_opacity

    elseif strcmp(rel_condition,'r')
        
        save Classifier_scores_labels_jitter.mat subject_trial_index stacked_SVM_avg_scores stacked_SVM_all_models_scores ...
            stacked_SVM_all_models_labels stacked_SVM_consensus_labels base_model_scores base_model_names ...
            imresize_num data_epoch modality relevant_location consensus_labels_* group* subject_trial_opacity
        
    end
    
elseif isequal(modality, 'MRI')
    
    save Classifier_scores_labels_jitter.mat subject_trial_index stacked_SVM_avg_scores stacked_SVM_all_models_scores ...
        stacked_SVM_all_models_labels stacked_SVM_consensus_labels base_model_scores base_model_names ...
        imresize_num data_epoch modality relevant_location subject_trial_opacity consensus_labels_* group*

end

%% Plot Stacked Model Confidence Score Histogram

%Setup figure
figure
hold on

%Define title
if strcmp(rel_condition, 'i')

    title(['Stacked SVM Model - ', modality, ' - ', irrelevant_location], 'Interpreter', 'none')

else strcmp(rel_condition, 'r')
    
    title(['Stacked SVM Model - ', modality, ' - ', relevant_location], 'Interpreter', 'none')

end

xlabel('Stacked SVM Confidence Scores')
ylabel('Number of Trials')
xlim([-6,6])

%Bin edges
bin_edges = [-6:0.2:6];

%Plot histogram of average confidence scores
histogram(stacked_SVM_avg_scores,bin_edges,'FaceColor','b')

%Plot reference lines
plot([0 0],[0 500], 'r')
plot([0.5 0.5],[0 500], '--k')
plot([-0.5 -0.5], [0 500], '--k')
plot([0.75 0.75],[0 500], '--k')
plot([-0.75 -0.75], [0 500], '--k')
plot([1 1],[0 500], '--k')
plot([-1 -1], [0 500], '--k')   

%Write text
text(2, 475, ['Perceived 0 = ',num2str(sum(stacked_SVM_avg_scores > 0))])
text(2, 450, ['Not perceived 0 = ',num2str(sum(stacked_SVM_avg_scores < 0))])

text(2, 425, ['Perceived 0.5 = ',num2str(sum(stacked_SVM_avg_scores > 0.5))])
text(2, 400, ['Not perceived -0.5 = ',num2str(sum(stacked_SVM_avg_scores < -0.5))])

text(-6, 475, ['Perceived 0.75 = ',num2str(sum(stacked_SVM_avg_scores > 0.75))])
text(-6, 450, ['Not perceived -0.75 = ',num2str(sum(stacked_SVM_avg_scores < -0.75))])

text(-6, 425, ['Perceived 1 = ',num2str(sum(stacked_SVM_avg_scores > 1))])
text(-6, 400, ['Not perceived -1 = ',num2str(sum(stacked_SVM_avg_scores < -1))])
    
%Save figure
cd(save_dir)
savefig(['Stacked_SVM_model_',modality,'_confidence_score_histogram.fig'])
close

%% Plot Eyelink Timecourses - Predicted Irrelevant vs True Class Relevant

%Load the report eyelink group data [time x subjects]
cd(report_data_dir)
load(fullfile(report_data_dir,report_data_filename));

%Exclude trials with blink at stimulus onset? Yes = 1; No = 0
exclude_blink_trials = 1;

%Run plotting function
eyelink_timecourse_predicted_irrelevant_jitter_function(save_dir,relevant_location, irrelevant_location, modality, group_CP_pupil_data,...
    group_CnP_pupil_data, group_CP_blink_data, group_CnP_blink_data, group_CP_microsac_data, group_CnP_microsac_data, noreport_data_pupil,...
    noreport_data_blink, noreport_data_microsac, exclude_blink_trials)
