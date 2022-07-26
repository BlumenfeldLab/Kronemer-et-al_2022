%% ERP average permutation stats topoplot - CP/CnP and PP/PnP T-values

%This code will take the average T-values for statistically significant
%clusters within specified time windows and generate a topoplot of those
%results. And adapated version of teh pop_topoplot function is used to only
%show the meaned time windows.

%Written by: Sharif Kronemer
%Date: 5/24/2021
%Modified: 7/12/2021

clear

%% Prompt

%Select run location
% prompt_1 = 'Running code local or server [l, s]: ';
% run_location = input(prompt_1,'s');
run_location = 's';

%Relevant trial type
prompt_2 = 'Data to analyze [CP,CnP,CP_minus_CnP]: ';
rel_data_type = input(prompt_2,'s');

%Irrelevant trial type
prompt_3 = 'Data to analyze [PP,PnP,PP_minus_PnP]: ';
irrel_data_type = input(prompt_3,'s');

%Permutation num
permutation_num = 5000;

%% Directories

%Data file name (Updated filename)
stat_folder_name = [rel_data_type,'_minus_',irrel_data_type,'_cent_quad_threshold_score_0.75_data_1501_3500ms_',num2str(permutation_num),'perm'];
%stat_folder_name = [irrel_data_type,'_minus_',rel_data_type,'_cent_quad_threshold_con_score_0.75_data_1501_3500ms_',num2str(permutation_num),'perm'];

%Save directory
if isequal(run_location, 's')
    
    addpath('/mnt/Data6/HNCT sEEG Study/sEEG Analysis/Analysis Code/eeglab14_0_0b')
    %addpath(genpath('/mnt/Data6/HNCT sEEG Study/sEEG Analysis/Analysis Code/eeglab14_0_0b/functions'))        
    addpath(genpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/EEG Analysis/eeglab14_0_0b/functions'))%/popfunc'))

    %Load EEGLab Template
    load('/mnt/Data6/HNCT sEEG Study/sEEG Analysis/Analysis Code/eeglab_prepared_blank.mat');
    
    %Save directory
    save_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/EEG Analysis/Permutation Analysis/Voltage/Cluster',...
        stat_folder_name);
    
    %Photogrammetry directory
    Photo_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/EEG Analysis/eeglab14_0_0b/sample_locs/GSN-HydroCel-257.sfp';    
    
elseif isequal(run_location, 'l')
    

end    

%% Find Stat Significant Clusters

%Plotting Method [t, voltage, constant]
plotting_method = 'voltage';

%Enter figure folder
cd(save_dir)
load(['permutation_stat_cluster_sumt_rand_005_',num2str(permutation_num),'_perm.mat'])

%If voltage plotting avg over subjects
if isequal(plotting_method,'voltage')
   
    %Average subject data
    group_avg_baselined_data = nanmean(group_pc_baselined_data,3);
    
end

%Create empty matrix channels by time
empty_matrix = zeros(size(group_pc_baselined_data,1),size(group_pc_baselined_data,2));

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
        if isequal(plotting_method, 'constant')
            
            %Give post clusters constant value - find the voxel all volume idx
            empty_matrix(sig_grid_voxels,time) = 0.7;
        
        %T-value
        elseif isequal(plotting_method, 't')
            
            %T-value
            empty_matrix(sig_grid_voxels,time) = t_orig(sig_grid_voxels,time);
            
        %Voltage
        elseif isequal(plotting_method, 'voltage')
            
            %Voltage
            empty_matrix(sig_grid_voxels,time) = group_avg_baselined_data(sig_grid_voxels,time);

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
        if isequal(plotting_method, 'constant')
            
            %Give post clusters constant value - find the voxel all volume idx
            empty_matrix(sig_grid_voxels,time) = -0.7;
        
        %T-value
        elseif isequal(plotting_method, 't')
            
            %T-value
            empty_matrix(sig_grid_voxels,time) = t_orig(sig_grid_voxels,time);
            
        %Voltage
        elseif isequal(plotting_method, 'voltage')
            
            %Voltage
            empty_matrix(sig_grid_voxels,time) = group_avg_baselined_data(sig_grid_voxels,time);            

        end
        
    end

end

%% Average T-values over ERP Windows 

%ERP time windows
N100 = [75:125];
VAN = [175:225];
P2_N2 = [275:325];
P3 = [350:650];

%ERP name
ERP_name = {'N100','VAN','P2_N2','P3'};

%Create cell of ERP time windows
ERP_times = {N100,VAN,P2_N2,P3};

%Define the center of the ERP time windows
ERP_center_times = [mean(N100),mean(VAN),mean(P2_N2),mean(P3)];

%Setup variable
all_ERP_matrix = [];

%Loop over ERPs
for type = 1:length(ERP_name)
    
    %Define current ERP name and time
    current_name = ERP_name{type};
    
    %Note: when considering the perm test on the -500 to 1500ms data, using
    %501 as the onset time
    current_time = ERP_times{type}+501; %Update the event times to correspond with epoch length
    
    %ERP mean t-value
    ERP_mean = nanmean(empty_matrix(:,current_time),2);
    
    %Add to all ERP matrix
    all_ERP_matrix = cat(2, all_ERP_matrix, ERP_mean);
    
end

%% Plot Topoplot of ERPs

%Set renderer to painter
set(0, 'DefaultFigureRenderer', 'painters');

%Define data
topoplot_data = all_ERP_matrix;

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

%Define the nose direction
eeglab_struct_voltage.chaninfo.nosedir = '+X';

%Define the data reference type (average reference)
eeglab_structure.ref = 'averef';

%Define the sampling rate
eeglab_struct_voltage.srate = 1000;

%Color map limits
if isequal(plotting_method,'t')
    
    color_bar_limits = [-4, 4];

elseif isequal(plotting_method,'voltage')
    
    color_bar_limits = [-2, 2];
    
end

%Array of time points to plot
time_points = ERP_center_times;

%Plot topoplot
pop_topoplot_ERP_timewindow(eeglab_struct_voltage, 1, time_points, [data_type,' voltage ERPs - N = ',num2str(size(group_pc_baselined_data,3))], [1 length(time_points)],0,... 
   'electrodes', 'on', 'style', 'map', 'intrad', [0.66], 'maplimits',color_bar_limts);

%Save figure
cd(save_dir)
savefig([rel_data_type,'_minus_',irrel_data_type,'_ERP_voltage_stat_',plotting_method,'_values_timepoint_topoplot.fig'])

close