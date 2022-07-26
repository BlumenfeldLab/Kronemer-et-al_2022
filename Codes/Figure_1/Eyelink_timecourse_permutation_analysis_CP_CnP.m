%% Eyelink Timecourse Permutation Analysis - CP/CnP

%This code will load the specified Eyelink data and run a temporal
%permutation test on the timecourse x subject matrix

%Written by: Sharif I. Kronemer
%Date: 7/4/2021

clear

%% Run Location

%Select run location
prompt_1 = 'Running code local or server [l, s]: ';
run_location = input(prompt_1,'s');

%Location relevant
prompt_3 = 'Relevant location condition [c, q, combine]: ';
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

%With trials with blink rejected
prompt_6 = 'Blink trials rejected [y,n]: ';
blink_trials_rejected = input(prompt_6,'s');

%Stimulus opacity
opacity = 'threshold';

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
            ['CP_CnP_',rel_save_name,'_',opacity,'_',modality_name]);
    
    else
        
        save_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/EyeLink Analysis/Permutation Analysis',...
            ['CP_CnP_',rel_save_name,'_',opacity,'_',modality_name,'_no_blink_rejections']);   

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
    
%% Load group Eyelink data
%Note: this subject data matrix is created with the
%Create_EyeLink_Group_CP_CnP_NRP_Dataset.m script

cd(data_dir)

if isequal(blink_trials_rejected,'y')
    
    load(['Group_eyelink_CP_CnP_',rel_save_name,'_',opacity,'_',modality_name,'_data.mat'])

else 
    
    load(['Group_eyelink_CP_CnP_',rel_save_name,'_',opacity,'_',modality_name,'_no_blink_rejections_data.mat'])

end

%% Run Permutation Analysis

%Data type list
data_list = {'pupil','blink','microsac'};

%Loop over data types (Pupil, blink, and microsaccade)
for data = 1:length(data_list)
    
    %Define current data type
    data_type = data_list{data};
    
    disp(['**Running - ',data_type])
    
    %Rename data to generic variable
    CP_eyelink_data = eval(['group_CP_',data_type,'_data']);
    CnP_eyelink_data = eval(['group_CnP_',data_type,'_data']);
    
    %% Main and Subtraction Epoch Subtraction and Baseline 

    disp('Baseline data')

    %Specify baseline period 
    baseline_window = 5001:6000;

    %CP Trials

    %Calculate baseline values [subjects]
    group_CP_baseline = squeeze(nanmean(CP_eyelink_data(baseline_window,:),1));

    %Convert baseline [time x subjects]
    group_CP_baseline = repmat(group_CP_baseline,[12000,1]);

    %Subtract baseline from main data
    group_CP_baselined_data = CP_eyelink_data - group_CP_baseline;

    %CnP Trials

    %Calculate baseline values [subjects]
    group_CnP_baseline = squeeze(nanmean(CnP_eyelink_data(baseline_window,:),1));

    %Convert baseline [time x subjects]
    group_CnP_baseline = repmat(group_CnP_baseline,[12000,1]);

    %Subtract baseline from main data
    group_CnP_baselined_data = CnP_eyelink_data - group_CnP_baseline;

    %% Run Permutation Tests

    disp('Running permutation test')

    %Restore zeros at stimulus presentaiton - Replace stim time values to
    %0; Note: original load data excludes trials with blink at stim time,
    %therefore the blink vectors should be 0 6001-6050ms. but when you
    %baseline data the 0s become non-zero values. The script below restores
    %those 0s for only the blink vector. 
    if isequal(data_type,'blink')
        
        group_CP_baselined_data(6001:6050,:) = 0;
        group_CnP_baselined_data(6001:6050,:) = 0;
        
    end
    
    
    tic

    %CP vs CnP Testing
    [clusters, pval, t_sums, permutation_distribution] = permutest_TimeCourses(group_CP_baselined_data, group_CnP_baselined_data, dependent_samples, ...
        p_threshold, num_permutations, two_sided);

    toc 

    %Find significant clusters pvalue < 0.05
    sig_clust = find(pval < 0.05);

    %Find the significant time points
    sig_time_pts = sort([clusters{sig_clust}]);

    %Save output
    cd(save_dir)
    save([data_type,'_timecourse_cluster_CP_vs_CnP.mat'],'clusters','pval','t_sums','permutation_distribution','sig_clust','sig_time_pts');

end
