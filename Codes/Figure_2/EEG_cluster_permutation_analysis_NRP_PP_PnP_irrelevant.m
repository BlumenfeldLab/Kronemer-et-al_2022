%% EEG Electrode Based Permutation Analysis - PP and PnP Task Irrelevant No report paradigm

%This code will:
%(1) Load EEG data
%(2) Baseline data
%(3) Cluster based permutation analysis 
%(4) Generate topoplots of significant clusters

%Written by: Sharif I. Kronemer
%Date: 2/24/2021
%Modified: 6/9/2021

clear

%% Run Location

%Select run location
%prompt_1 = 'Running code local or server [l, s]: ';
%run_location = input(prompt_1,'s');
run_location = 's';

%Trial type
prompt_2 = 'Data to analyze [PP,PnP,PP_minus_PnP]: ';
data_type = input(prompt_2,'s');

%Confidence score
prompt_3 = 'Confidence score [0,0.25,0.5,0.75]: ';
confidence_score = input(prompt_3,'s');

%% Parameters

%Number of permutations
num_permutation = 5000;

%Center and quadrant name variable (e.g., 'cent','quad', and 'cent_quad')
location_name = 'cent_quad';

%Stimulus opacity type (threshold/blank)
stim_opacity = 'blank';

%Irrelevant/Relevant conditions (relevant/irrelevant/rel_irrel)
task_condition = 'rel_irrel';

%% Directories and Variable Names

%Variable name
folder_name = [location_name,'_', stim_opacity];

if isequal(run_location, 's')
    
    %Add behavioral analysis path
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/Behavioral Analysis')
    
    %Add SPM and supplementary functions to path
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Permutation Statistical Analysis/Permutation functions')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/EEG Analysis/Permutation Statistical Analysis')
    addpath('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Time courses/Supplementary functions/fMRI Analysis/Code for Percent Change Maps/subRoutines')

    addpath('/mnt/Data6/HNCT sEEG Study/sEEG Analysis/Analysis Code/eeglab14_0_0b')
    addpath(genpath('/mnt/Data6/HNCT sEEG Study/sEEG Analysis/Analysis Code/eeglab14_0_0b/functions'))        

    %Load EEGLab Template
    load('/mnt/Data6/HNCT sEEG Study/sEEG Analysis/Analysis Code/eeglab_prepared_blank.mat');

    %Photogrammetry directory
    Photo_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/EEG Analysis/eeglab14_0_0b/sample_locs/GSN-HydroCel-257.sfp';
    
    %Data directory
    data_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Data/EEG/Voltage';
    
    %Save directory
    save_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/EEG Analysis/Permutation Analysis/Voltage/Cluster',...
        [data_type,'_',folder_name,'_con_score_',num2str(confidence_score),'_data_1501_3500ms_',num2str(num_permutation),'perm']);
    
    %Neighborhood matrix directory
    hood_dir = '/mnt/Data6/HNCT sEEG Study/Testing Documents';
    
elseif isequal(run_location, 'l')
     
end

%Make save directory
mkdir(save_dir)

%% Load group voltage data

cd(data_dir)
%load(['Group_EEG_voltage_irrelevant_PP_PnP_',folder_name,'_score_thres_',num2str(confidence_score),'_data.mat'])
load(['Group_EEG_voltage_',task_condition,'_PP_PnP_',folder_name,'_score_thres_',num2str(confidence_score),'_data.mat'])

%% Subtraction and Baseline 

disp('Baselining data from prestim period')

%Specify baseline period 
baseline_window = 1001:2000; 

%Subtraction
if isequal(data_type, 'PP_minus_PnP')
    
    %Check if the number of subjects is equal for PP and PnP group
    %variables
    if size(PP_epochs_subjects_list,1)> size(PnP_epochs_subjects_list,1) %If PP is greater
        
        %Delete subjects not found in the PnP class
        group_PP_EEG_voltage_data(:,:,not(ismember(PP_epochs_subjects_list,PnP_epochs_subjects_list))) = [];
        
        group_PnP_EEG_voltage_data(:,:,[41,56]) = [];
        
    elseif size(PP_epochs_subjects_list,1) < size(PnP_epochs_subjects_list,1) %If PnP is greater
        
        %Delete subjects not found in the PP class
        group_PnP_EEG_voltage_data(:,:,not(ismember(PnP_epochs_subjects_list,PP_epochs_subjects_list))) = [];
        
    end
    
    %Subtract main data from subtract data
    group_pc_data = group_PP_EEG_voltage_data - group_PnP_EEG_voltage_data;

    %Calculate baseline values [voxel x time x subjects]
    group_pc_baseline = squeeze(nanmean(group_pc_data(:,baseline_window,:),2));
    
    %Convert baseline [channel x time x subjects]
    group_pc_baseline = permute(repmat(group_pc_baseline,[1,1,4001]),[1,3,2]);
    
    %Subtract baseline from main data
    group_pc_baselined_data = group_pc_data - group_pc_baseline;

%     %Loop over time points
%     for time = 1:size(group_CP_EEG_voltage_data,2)
% 
%         %Subtract baseline from all voxels at all times
%         group_pc_baselined_data(:,time,:) = squeeze(group_pc_data(:,time,:)) - group_pc_baseline;
% 
%     end

%No subtraction    
elseif isequal(data_type, 'PP')
    
    %Calculate baseline values [channel x subjects]
    group_pc_baseline = squeeze(nanmean(group_PP_EEG_voltage_data(:,baseline_window,:),2));
    
    %Convert baseline [channel x time x subjects]
    group_pc_baseline = permute(repmat(group_pc_baseline,[1,1,4001]),[1,3,2]);
    
    %Subtract baseline from main data
    group_pc_baselined_data = group_PP_EEG_voltage_data - group_pc_baseline;

%     %Loop over time points
%     for time = 1:size(group_CP_EEG_voltage_data,2)
% 
%         %Subtract baseline from all values
%         group_pc_baselined_data(:,time,:) = squeeze(group_CP_EEG_voltage_data(:,time,:)) - group_pc_baseline;
% 
%     end
    
elseif isequal(data_type, 'PnP')
    
    %Calculate baseline values [voxel x time x subjects]
    group_pc_baseline = squeeze(nanmean(group_PnP_EEG_voltage_data(:,baseline_window,:),2));

    %Convert baseline [channel x time x subjects]
    group_pc_baseline = permute(repmat(group_pc_baseline,[1,1,4001]),[1,3,2]);
    
    %Subtract baseline from main data
    group_pc_baselined_data = group_PnP_EEG_voltage_data - group_pc_baseline;
    
%     %Loop over time points
%     for time = 1:size(group_CnP_EEG_voltage_data,2)
% 
%         %Subtract baseline from all values
%         group_pc_baselined_data(:,time,:) = squeeze(group_CnP_EEG_voltage_data(:,time,:)) - group_pc_baseline;
% 
%     end
    
end

%% Manage Temporal Resolution

%Select the samples of time to complete the permutation test
test_window = 1501:3500;
group_pc_baselined_data = group_pc_baselined_data(:,test_window,:);

% %Bin time 
% group_pc_baselined_data_downsample = discretize(group_pc_baselined_data(1,:,1),100);
% 
% end

%% Neighborhood Matrix

%Load neighborhood matrix    
load(fullfile(hood_dir,'net_neighborhood_matrix_257.mat'))

%% Run Permutation Test

disp('Running Permutation Test')

tic

%INPUTs: 3D data matrix, 2D Neighborhood matrix
[pval, t_orig, clust_info, seed_state, est_alpha, mn_clust_mass] = EEG_clust_perm1_iceeg_sumt_rand(group_pc_baselined_data, net_neighborhood, num_permutation, 0.05, 0, 0.05, 2, [], 0);

toc

%Save output
cd(save_dir)
save(['permutation_stat_cluster_sumt_rand_005_',num2str(num_permutation),'_perm.mat'],'pval','t_orig','clust_info','seed_state','est_alpha','mn_clust_mass' ...
,'group_pc_baselined_data');

%% Plot Topoplot

%Plotting Method
plotting_method = 'v';%'t';

%Mean voltage over subjects
group_voltage_data = nanmean(group_pc_baselined_data,3);

%Enter figure folder
cd(save_dir)
load(['permutation_stat_cluster_sumt_rand_005_',num2str(num_permutation),'_perm.mat'])

%Create empty matrix channels by time
%empty_matrix = zeros(257,4001);
empty_matrix = zeros(257,2000);
%empty_matrix = zeros(257,500);

%Find significant clusters pvalue < 0.05 - Negative and positive clusters independently
pos_clust = find(clust_info.pos_clust_pval < 0.05);
neg_clust = find(clust_info.neg_clust_pval < 0.05);

%Loop over postive clusters
for clust = 1:length(pos_clust)
    
    %Loop over time
    for time = 1:size(empty_matrix,2)
              
        %Find postive significant cluster number/idx in voxel x time matrix
        [sig_grid_voxels, col] = find(clust_info.pos_clust_ids(:,time) == pos_clust(clust));
        
        %Constant
        if isequal(plotting_method, 'c')
            
            %Give post clusters constant value - find the voxel all volume idx
            empty_matrix(sig_grid_voxels,time) = 0.7;
        
        %T-value
        elseif isequal(plotting_method, 't')
            
            %T-value
            empty_matrix(sig_grid_voxels,time) = t_orig(sig_grid_voxels,time);
        
        %Voltage
        elseif isequal(plotting_method, 'v')
            
            %Sig voltage
            empty_matrix(sig_grid_voxels,time) = group_voltage_data(sig_grid_voxels,time);
            
        end
            
    end

end

%Loop over postive clusters
for clust = 1:length(neg_clust)
    
    %Loop over time
    for time = 1:size(empty_matrix,2)
              
        %Find postive significant cluster number/idx in voxel x time matrix
        [sig_grid_voxels, col] = find(clust_info.neg_clust_ids(:,time) == neg_clust(clust));
        
        %Constant
        if isequal(plotting_method, 'c')
            
            %Give post clusters constant value - find the voxel all volume idx
            empty_matrix(sig_grid_voxels,time) = -0.7;
        
        %T-value
        elseif isequal(plotting_method, 't')
            
            %T-value
            empty_matrix(sig_grid_voxels,time) = t_orig(sig_grid_voxels,time);

        %Voltage
        elseif isequal(plotting_method, 'v')
            
            %Sig voltage
            empty_matrix(sig_grid_voxels,time) = group_voltage_data(sig_grid_voxels,time);
        
        end
        
    end

end

%% Plot Voltage Topoplots

%Topoplot timing cells
start_time_cell = [-100];%[-50,-100,-250];
end_time_cell = [1499];%[1150, 250, 0];
increment_cell = [100];%[50, 25, 25];

%Loop over timing options
for timing = 1:length(start_time_cell)

    %Define time points for topoplots
    start_time = start_time_cell(timing);
    end_time = end_time_cell(timing);
    increment = increment_cell(timing);    
    
    %Define data
    topoplot_data = empty_matrix;

    %Setup structure for topoplot
    eeglab_struct_voltage = eeglab_prepared_blank;

    %Set channel locations
    eeglab_struct_voltage.chanlocs = readlocs(Photo_dir); 
    eeglab_struct_voltage.urchanlocs = readlocs(Photo_dir);

    %Add data to structure
    eeglab_struct_voltage.data = topoplot_data;

    %Define the number of trials or subjects (EEGLAB averages over the 3rd dimension)
    eeglab_struct_voltage.trials = [];

    %Define the number of time points
    eeglab_struct_voltage.pnts = size(topoplot_data,2);

    %Define the number of channels
    eeglab_struct_voltage.nbchan = size(topoplot_data,1);

    %Define epoch onset and offset time (s)
    eeglab_struct_voltage.xmin = -0.5;%-2; 

    %Define epoch onset and offset time
    eeglab_struct_voltage.xmax = 1.5;%2; 
    
    %Define stimulus onset time (sample number)
    eeglab_struct_voltage.stim_onset = 501;

    %Define the nose direction
    eeglab_struct_voltage.chaninfo.nosedir = '+X';

    %Define the data reference type (average reference)
    eeglab_structure.ref = 'averef';

    %Define the sampling rate
    eeglab_struct_voltage.srate = 1000;

    %Color map limits
    if isequal(plotting_method, 'v')
        
        color_bar_limts = [-2, 2];
        
        data_name = 'voltage';

    else
        
        color_bar_limts = [-5, 5];
        
        data_name = 't-values';
        
    end

    %Array of time points to plot
    time_points = [start_time:increment:end_time];

    %Plot topoplot
    %pop_topoplot_voltage(eeglab_struct_voltage, 1, time_points, [data_type,' voltage ERPs - N = ',num2str(size(group_pc_baselined_data,3))], [1 length(time_points)],0,... 
    %    'electrodes', 'on', 'style', 'map', 'intrad', [0.7], 'maplimits',color_bar_limts);
    pop_topoplot_voltage(eeglab_struct_voltage, 1, time_points, [data_type,' voltage ',data_name, ' ERPs - N = ',num2str(size(group_pc_baselined_data,3))], [1 length(time_points)],0,... 
        'electrodes', 'on', 'style', 'map','gridscale',50, 'emarker',{'.','k',3},'plotchans',[], 'intrad', [0.5], 'maplimits',color_bar_limts);
    
    %Save figure
    cd(save_dir)
    if isequal(plotting_method,'v')
        
        savefig([data_type, '_voltage_topoplot_',num2str(start_time),'_',num2str(end_time),'ms.fig'])

    else
        
        savefig([data_type, '_voltage_t_values_topoplot_',num2str(start_time),'_',num2str(end_time),'ms.fig'])
    
    end
    
    close

end