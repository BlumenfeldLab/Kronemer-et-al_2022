%% Eyelink Timecourse Permutation Analysis

%This code will load the specified Eyelink data and run a temporal
%permutation test on the timecourse x subject matrix

%Written by: Sharif I. Kronemer
%Date: 3/10/2021
%Modified: 5/15/2022

clear

%% Run Location

%Select run location
prompt_1 = 'Running code local or server [l, s]: ';
run_location = input(prompt_1,'s');

%Select relevant or irrelevant condition
prompt_6 = 'Relevant or irrelvant condition [relevant, irrelevant, rel_irrel]: ';
relevant_condition = input(prompt_6,'s');

%Location relevant
prompt_3 = 'No-report relevant location [c, q, combine]: ';
relevant_location = input(prompt_3, 's');

%Center Relevant
if isequal(relevant_location, 'c')
    
    relevant_location = 'Center Relevant';
    rel_save_name = 'cent';

%Quadrant Relevant
elseif isequal(relevant_location, 'q')
    
    relevant_location = 'Quadrant Relevant';
    rel_save_name = 'quad';
    
%Combine Center and Quadrant Relevant   
elseif isequal(relevant_location, 'combine')
    
    relevant_location = 'Center and Quadrant Relevant';
    rel_save_name = 'cent_quad';
    
end

%Imaging modality
prompt_5 = 'Modality type [EEG, MRI, both]: ';
modality = input(prompt_5,'s');

%EEG
if isequal(modality,'EEG')
   
    modality_name = 'EEG';
    
%MRI
elseif isequal(modality,'MRI')
    
    modality_name = 'MRI';
    
%EEG and MRI
elseif isequal(modality,'both')
    
    modality_name = 'EEG_MRI';
    
end

%Confidence threshold value
prompt_4 = 'Confidence threshold value [0,0.25,0.5,0.75]: ';
confidence_score = input(prompt_4,'s');

%With trials with blink rejected
prompt_6 = 'Blink trials rejected [y,n]: ';
blink_trials_rejected = input(prompt_6,'s');

%Stimulus opacity
prompt_7 = 'Stimulus opacity [threshold, blank]: ';
opacity = input(prompt_7,'s');

%% Directories 

if isequal(run_location, 's')
    
    %Add behavioral analysis path
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/Behavioral Analysis')
    
    %Add paths
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Permutation Statistical Analysis/Permutation functions')
    addpath('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Time courses/Supplementary functions/fMRI Analysis/Code for Percent Change Maps/subRoutines')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Permutation Statistical Analysis')

    %Data directory
    data_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Data/EyeLink';
    
    %Save directory
    if isequal(blink_trials_rejected,'y')
        
        save_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/EyeLink Analysis/Permutation Analysis',...
            ['PP_PnP_',relevant_condition,'_',rel_save_name,'_',opacity,'_',modality_name,'_confid_score_',num2str(confidence_score)]);

    else
        
        save_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/EyeLink Analysis/Permutation Analysis',...
            ['PP_PnP_',relevant_condition,'_',rel_save_name,'_',opacity,'_',modality_name,'_confid_score_',num2str(confidence_score),'_no_blink_rejections']);
        
    end   
        
elseif isequal(run_location, 'l')
       
end

%Make save directory
mkdir(save_dir)

%% Parameters

%Number of permutations
num_permutations = 5000;

%Are samples dependent (default is true)
dependent_samples = 'true';

%Define alpha threshold
p_threshold = 0.05;

%Two-sided
two_sided = 'true';
    
%% Load group voltage data

cd(data_dir)

if isequal(blink_trials_rejected,'y')

    load(['Group_eyelink_',relevant_condition,'_PP_PnP_',rel_save_name,'_',opacity,'_',modality_name,'_score_thres_',num2str(confidence_score),'_data.mat'])

else
    
    load(['Group_eyelink_',relevant_condition,'_PP_PnP_',rel_save_name,'_',opacity,'_',modality_name,'_score_thres_',num2str(confidence_score),'_no_blink_rejections_data.mat'])

end

%% Run Permutation Analysis

%Data type list
data_list = {'pupil','blink','microsac'};

%Loop over data types
for data = 1:length(data_list)
    
    %Define current data type
    data_type = data_list{data};
    
    disp(['**Running - ',data_type])
    
    %Rename data to generic variable
    PP_eyelink_data = eval(['group_PP_',data_type,'_data']);
    PnP_eyelink_data = eval(['group_PnP_',data_type,'_data']);
    
    %Check if PP and PnP are the same size
    if size(PP_eyelink_data,2) > size(PnP_eyelink_data,2)
        
       PP_eyelink_data(:,~ismember(PP_epochs_subjects_list,PnP_epochs_subjects_list)) = [];
        
    elseif size(PP_eyelink_data,2) < size(PnP_eyelink_data,2)
        
       PnP_eyelink_data(:,~ismember(PnP_epochs_subjects_list,PP_epochs_subjects_list)) = [];
            
    end
    
    %% Main and Subtraction Epoch Subtraction and Baseline 

    disp('Baseline data')

    %Specify baseline period 
    baseline_window = 5001:6000;

    %CP Trials

    %Calculate baseline values [subjects]
    group_PP_baseline = squeeze(nanmean(PP_eyelink_data(baseline_window,:),1));

    %Convert baseline [time x subjects]
    group_PP_baseline = repmat(group_PP_baseline,[12000,1]);

    %Subtract baseline from main data
    group_PP_baselined_data = PP_eyelink_data - group_PP_baseline;

    %CnP Trials

    %Calculate baseline values [subjects]
    group_PnP_baseline = squeeze(nanmean(PnP_eyelink_data(baseline_window,:),1));

    %Convert baseline [time x subjects]
    group_PnP_baseline = repmat(group_PnP_baseline,[12000,1]);

    %Subtract baseline from main data
    group_PnP_baselined_data = PnP_eyelink_data - group_PnP_baseline;

    %% Run Permutation Tests

    disp('Running permutation test')
    
    %Restore zeros at stimulus presentaiton - Replace stim time values to
    %0; Note: original load data excludes trials with blink at stim time,
    %therefore the blink vectors should be 0 6001-6050ms. but when you
    %baseline data the 0s become non-zero values. The script below restores
    %those 0s for only the blink vector. 
    if isequal(data_type,'blink')
        
        group_PP_baselined_data(6001:6050,:) = 0;
        group_PnP_baselined_data(6001:6050,:) = 0;
        
    end

    tic

    %CP vs CnP Testing
    [clusters, pval, t_sums, permutation_distribution] = permutest_TimeCourses(group_PP_baselined_data, group_PnP_baselined_data, dependent_samples, ...
        p_threshold, num_permutations, two_sided);

    toc 

    %Find significant clusters pvalue < 0.05
    sig_clust = find(pval < 0.05);

    %Find the significant time points
    sig_time_pts = sort([clusters{sig_clust}]);

    %Save output
    cd(save_dir)
    save([data_type,'_timecourse_cluster_PP_vs_PnP.mat'],'clusters','pval','t_sums','permutation_distribution','sig_clust','sig_time_pts');

end