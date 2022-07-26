%% ERP average permutation stats topoplot - PP/PnP T-values

%This code will take the average T-values for statistically significant
%clusters within specified time windows and generate a topoplot of those
%results. And adapated version of teh pop_topoplot function is used to only
%show the meaned time windows.

%Written by: Sharif Kronemer
%Date: 5/24/2021
%Modified: 7/2/2022

clear

%% Prompt

%Select run location
% prompt_1 = 'Running code local or server [l, s]: ';
% run_location = input(prompt_1,'s');
run_location = 's';

%Plotting Method [Constant value = c, t-value = t; significant times voltage = v_sig; all times voltage = v_all]
prompt_4 = 'Plotting method [c,t,v_sig,v_all]: ';
plotting_method = input(prompt_4,'s');

%% Directories
    
addpath('/mnt/Data6/HNCT sEEG Study/sEEG Analysis/Analysis Code/eeglab14_0_0b')
%addpath(genpath('/mnt/Data6/HNCT sEEG Study/sEEG Analysis/Analysis Code/eeglab14_0_0b/functions'))        
addpath(genpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/EEG Analysis/eeglab14_0_0b/functions'))%/popfunc'))

%Load EEGLab Template
load('/mnt/Data6/HNCT sEEG Study/sEEG Analysis/Analysis Code/eeglab_prepared_blank.mat');

%Data directory
data_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Data/EEG/Voltage');

%Photogrammetry directory
Photo_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/EEG Analysis/eeglab14_0_0b/sample_locs/GSN-HydroCel-257.sfp';    

%% Run Analysis

%Data cell
data_types_cell = {'PP','PnP','PP_minus_PnP'};

%Loop over data types
for type = 1:length(data_types_cell)
    
    %Select data type
    data_type = data_types_cell{type};
    
    disp(['***Running ',data_type])

    %Data file name (Update filename below)
    stat_folder_name = [data_type,'_cent_quad_threshold_con_score_0.75_data_1501_3500ms_5000perm'];
    %stat_folder_name = [data_type,'_cent_quad_blank_con_score_0.75_data_1501_3500ms_100perm'];
    %stat_folder_name = [data_type,'_microsaccade_cent_quad_blank_data_1501_3500ms_100perm'];
    %stat_folder_name = [data_type,'_cent_quad_jitter_matched_con_score_0.75_data_1501_3500ms_5000perm'];

    %Save directory
    save_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/EEG Analysis/Permutation Analysis/Voltage/Cluster',...
        stat_folder_name);

    %% Average Voltage over ERP Windows 

    %Data file name (Update filename below)
    data_file_name = 'Group_EEG_voltage_irrelevant_PP_PnP_cent_quad_threshold_score_thres_0.75_data.mat';
    %data_file_name = 'Group_EEG_voltage_rel_irrel_PP_PnP_cent_quad_blank_score_thres_0.75_data.mat';
    %data_file_name = 'Group_EEG_voltage_rel_irrel_microsaccade_cent_quad_blank_data.mat';
    %data_file_name = 'Group_EEG_voltage_rel_irrel_PP_PnP_cent_quad_jitter_matched_score_thres_0.75_data.mat';

    %Load group voltage data
    load(fullfile(data_dir,data_file_name))

    %Define/Setup data type
    if isequal(data_type,'PP')

        group_plotting_data = group_PP_EEG_voltage_data;

    elseif isequal(data_type,'PnP')

        group_plotting_data = group_PnP_EEG_voltage_data;

    elseif isequal(data_type,'PP_minus_PnP')

        group_plotting_data = nanmean(group_PP_EEG_voltage_data,3) - nanmean(group_PnP_EEG_voltage_data,3);

    end

    %Count subjects
    num_subjects = size(group_plotting_data,3);

    %Average over subjects 
    group_plotting_data = nanmean(group_plotting_data,3);

    %Cut time dimesion to the stat test window
    group_plotting_data = group_plotting_data(:,1501:3500);

    %% Find Stat Significant Clusters

    %Enter figure folder
    cd(save_dir)
    load('permutation_stat_cluster_sumt_rand_005_5000_perm.mat')

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
            if isequal(plotting_method, 'c')

                %Give post clusters constant value - find the voxel all volume idx
                empty_matrix(sig_grid_voxels,time) = 1;

            %T-value
            elseif isequal(plotting_method, 't')

                %T-value
                empty_matrix(sig_grid_voxels,time) = t_orig(sig_grid_voxels,time);

            %Voltage
            elseif isequal(plotting_method, 'v_all') || isequal(plotting_method, 'v_sig')

                empty_matrix(sig_grid_voxels,time) = group_plotting_data(sig_grid_voxels,time);


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
                empty_matrix(sig_grid_voxels,time) = -1;

            %T-value
            elseif isequal(plotting_method, 't')

                %T-value
                empty_matrix(sig_grid_voxels,time) = t_orig(sig_grid_voxels,time);

            %Voltage
            elseif isequal(plotting_method, 'v_all') || isequal(plotting_method, 'v_sig')

                empty_matrix(sig_grid_voxels,time) = group_plotting_data(sig_grid_voxels,time); 

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

        %ERP mean values within time window
        if isequal(plotting_method,'v_all')
        
            ERP_mean = nanmean(group_plotting_data(:,current_time),2);
        
        else
            
            ERP_mean = nanmean(empty_matrix(:,current_time),2);
    
        end
    
        %Find channels with sig times
        ERP_channels = empty_matrix(:, current_time);
        ERP_channels(ERP_channels ~= 0) = 1; %Replace all non-zeros with 1
        ERP_channels = sum(ERP_channels,2);
        ERP_channels = find(ERP_channels > round(length(current_time).*0.5));
    
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

        color_bar_limts = [-2, 2]; %[-8, 8];

    elseif isequal(plotting_method,'c')

        color_bar_limts = [-1, 1];

    elseif isequal(plotting_method,'v_all') || isequal(plotting_method,'v_sig')

        color_bar_limts = [-2, 2];

    end

    %Array of time points to plot
    time_points = ERP_center_times;

    %Plot topoplot
    %pop_topoplot_ERP_timewindow(eeglab_struct_voltage, 1, time_points, [data_type,' voltage ERPs - N = ',num2str(size(group_pc_baselined_data,3))], [1 length(time_points)],0,... 
    %    'electrodes', 'on', 'style', 'map', 'intrad', [0.66], 'maplimits',color_bar_limts);    
    %pop_topoplot_ERP_timewindow(eeglab_struct_voltage, 1, time_points, [data_type,' voltage stat ERPs - N = ',num2str(size(group_pc_baselined_data,3))], [1 length(time_points)],0,... 
    %    'electrodes', 'on', 'style', 'map','gridscale',50, 'emarker',{'.','k',3},'plotchans',[],'intrad', [0.5], 'maplimits',color_bar_limts);
    %pop_topoplot_ERP_timewindow(eeglab_struct_voltage, 1, time_points, [data_type,' voltage stat ERPs - N = ',num2str(size(group_pc_baselined_data,3))], [1 length(time_points)],0,... 
    %    'electrodes', 'on', 'style', 'map','gridscale',50, 'emarker',{'.','k',3},'emarker2',{ERP_channels,'o','w',2},'plotchans',[],'intrad', [0.5], 'maplimits',color_bar_limts);
    pop_topoplot_ERP_timewindow(eeglab_struct_voltage, 1, time_points, [data_type,' voltage stat ERPs - N = ',num2str(size(group_pc_baselined_data,3))], [1 length(time_points)],0,... 
        'electrodes', 'on', 'style', 'map','gridscale',50, 'emarker',{'.','k',3},'plotchans',[],'intrad', [0.5],'maplimits',color_bar_limts);

    %Save figure
    cd(save_dir)

    if isequal(plotting_method,'t')

        savefig([data_type,'_ERP_stat_sig_t_values_timepoint_topoplot.fig'])

    elseif isequal(plotting_method,'c')

        savefig([data_type,'_ERP_stat_sig_constant_timepoint_topoplot.fig'])

    elseif isequal(plotting_method,'v_sig')

        savefig([data_type,'_ERP_stat_sig_voltage_timepoint_topoplot.fig'])
        
    elseif isequal(plotting_method,'v_all')

        savefig([data_type,'_ERP_stat_all_voltage_timepoint_topoplot.fig'])
      
    end

    close

end
