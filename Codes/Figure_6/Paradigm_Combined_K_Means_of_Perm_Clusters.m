%% Combined Report and No Report Dataset K-means clustering of permutation clusters

%This code will take the permutation statistical results and apply for
%clustering in space by the correlation in percent change BOLD timecourses
%within the post-stimulus period (pre questions). Importantly, the code
%will only consider voxels in the k-means clustering that was found
%statistically significant within the same post-stim timeframe, regardless
%to whether that voxel was significant for only 1 or all time points in
%the timeframe. 

% Name: Sharif Kronemer
% Date: 07/2/2020
% Modified: 6/24/2021

% uses pearon's correlation to cluster the voxel data.
% implemented using the kmeans matlab function https://www.mathworks.com/help/stats/kmeans.html

% Need matlab stats toolbox, marsbar (for creating the ROI)

% Input: percent change maps, clustering parameters
% Outputs:
%       -In the output directory:
%           -[Cluster_num_...] folder (based on number of clusters we
%           clustered with) which contains: 
%               -Clustering info: clusterinfo.mat
%               -Cluster ROIs
%               -(Cluster timecourses)
%           -Silhouette plots

% clusterinfo.mat: 
% -Clusters contain cluster indices (1 through cluster_num)
% -C contains the centroid locations (cluster_num by t)
% -sumd contains the within-cluster sums to centroid distances (cluster_num
% by 1)
% -D contains the distances from each point to each centroid (voxel_num by
% cluter_num)

clear

%% Prompts

%Run location
prompt1 = 'Run location (s/l): ';
run_location = input(prompt1, 's');

%Select Positive or Negative Clusters
prompt2 = 'Positive Clusters, Negative Clusters, or Both (p,n,b): ';
sig_clusters = input(prompt2, 's');

%% Parameters 

%Variable used for establishing directories and folder names
ROI_name = 'Whole Brain';

%% Add Functions to Path

%Local
if isequal(run_location, 'l')
    
%Server
elseif isequal(run_location, 's')
    
    %Add SPM and supplementary functions to path
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Permutation Statistical Analysis/Permutation functions')
    addpath('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Time courses/Supplementary functions/fMRI Analysis/Code for Percent Change Maps/subRoutines')
    addpath('/mnt/Data5/Sharif Data F/HNCT fMRI Study/spm12_radiological')
    addpath(genpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Extract voxel data/Supplementary functions'));
    addpath(['/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/marsbar-0.44']);
   
    %MNI brain
    %template_image = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/spm12_radiological/canonical/single_subj_T1.nii';
    template_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Subject Analysis/191AM/MRI Analysis/Extracted voxel data/Run 1 Movie/axialSlices'; 

    %Save dir
    save_dir = '/mnt/Data8/HNCT NRP and Report Only Paradigm Analyses/Analysis/MRI Analysis/K means Clusters';
        
    %Cluster directory
    stat_dir = '/mnt/Data8/HNCT NRP and Report Only Paradigm Analyses/Analysis/MRI Analysis/Percent Change Permutation Stats/Whole Brain CP minus CnP Avg 5000';
    
    %Voxel map
    %voxel_dir = ['/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Permutation Statistical Analysis/Brain Voxel Subset/', ROI_name];
    
    %MNI brain
    %ROI_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Subject Analysis/191AM/MRI Analysis/Extracted voxel data/Run 1 Movie/axialSlices'; 
    ROI_dir = ['/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Permutation Statistical Analysis/Brain Voxel Subset/', ROI_name]; %Cortical/Subcortical voxels
    
    %Data dir
    data_dir = '/mnt/Data8/HNCT NRP and Report Only Paradigm Analyses/Group Data/MRI/BOLD PC';
    
    %Cortical and subcortical ROI (Custom made cortical and subcortical ROI from MarsBaR and harvard atlas)
    subcortical_ROI_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/Subcortical ROI';
    cortical_ROI_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/Marsbar ROIs/Cortical ROI';

end

%% Parameters

%Clustering parameters
max_clusters = 10;  %max number of clusters, these many clusters will be attempted
replicate_num = 10; %will re-run clustering with new initial values
maxiter = 400; %number of iterations within clustering events (i.e., find new cluster centroid, recluster values, find new centroid,...etc.)
timeframe = 21:31; %the epoch in the precent change map data that will be clustered (e.g., 21s to 36s represents the time from face onset to the end of the 15s post-stim window)

%% Generate Group Percent Change Variable [time x voxel x subjects]

disp('Load BOLD PC Data')

cd(data_dir)
load('Group_PC_BOLD_CP_CnP_threshold_data.mat')

%Subtract main data from subtract data
avg_group_pc_data = nanmean(group_CP_BOLD_PC_data,3) - nanmean(group_CnP_BOLD_PC_data,3);
        
%clearvars group_CP_BOLD_PC_data group_CnP_BOLD_PC_data

%% Load Permutation Clusters and Select Significant Clusters (p < 0.05)

%Voxel Method defines if all or a subset of voxels will be considered in clustering
%0 = use all voxels 
%1 = use grid/ROI subset voxels
voxel_method = 1;

%Load permutation cluster information
cd(stat_dir)
load('pc_cluster_sumt_rand_005.mat')

%Find significant clusters pvalue < 0.05 - Negative and positive clusters independently
pos_clust = find(clust_info.pos_clust_pval < 0.05);
neg_clust = find(clust_info.neg_clust_pval < 0.05);

%Load ROI voxelInfo
cd(ROI_dir)
load(fullfile(ROI_dir,'voxelInfo.mat'));

% %Define XYZ position and select grid subset
% try
%     
%     %Voxel info
%     vXYZ = xyz.vXYZ;
%     
%     %Image matrix
%     %V = xyz.mat;
%     
% catch
%     
%     %Voxel info
%     vXYZ = voxelInfo;
%     
%     %Image matrix
%     %V = transform_matrix;
%     
% end

%Load MNI Brain
cd(template_dir)
load(fullfile(template_dir,'voxelInfo.mat'));

%Voxel info
temp_vXYZ = xyz.vXYZ; %MNI

%Image matrix
V = xyz.mat;

%Use grid voxels
if isequal(voxel_method, 1)
    
    disp('Cutting voxels to grid/ROI subset!')
    
    %Load voxel grid and adjacency matrix
    cd(ROI_dir)
    
    try
        
        %Whole brain grid
        load('brain_grid_voxel_info.mat')
    
    catch
        
        %ROI grid 
        load('ROI_grid_voxel_info.mat')
               
    end
    
    %Cut voxels positions to only those in the brain_grid_voxels matrix 
    vXYZ = vXYZ(:, brain_grid_voxels); %ROI

%Use all voxels
else
   
    disp('Using all voxels!')
    
end   

%Specify image size (should be 91x109x91);
sizeMat = [91,109,91]; 

%Get the list of index of useful voxels sizeMat(2:4), is the range of X, Y, Z indices. sizeMat(1) is time
grid_outInd = mat2ind(vXYZ,sizeMat); %ROI

%outInd = mat2ind(vXYZ,sizeMat);
outInd = mat2ind(temp_vXYZ,sizeMat); %MNI Brain

%% Average BOLD From Adjacent Voxels

%This for loop cycles through all time points per participant
for timepoint = 1:size(avg_group_pc_data,1) 

    %Reshape session voxel data to 91, 109, 91 for each time point and subject
    img = reshape(avg_group_pc_data(timepoint,:),[91, 109, 91]);

    %Loop over voxels
    for voxel = 1:size(vXYZ,2)
        
        %Use grid voxels and average over adjacent voxels
        if isequal(voxel_method,1)
            
            %Find the mean BOLD signal for all voxels within spatial
            %cluster, includes center and adjacent voxels - nanmean([adjacent voxels + center voxel])
            clustering_matrix(voxel,timepoint) = nanmean([img(all_adjacent_voxels_cell{voxel}); img(grid_outInd(voxel))]);
        
        %Use all voxels, no adjacent voxels to average
        else
            
            clustering_matrix(voxel,timepoint) = nanmean(img(grid_outInd(voxel)));
  
        end
        
    end

end

%% Select Subset of Significant Voxels

%Matrix of positive and negative clusters
all_sig_pos_voxels = [];
all_sig_neg_voxels = [];

%Loop over postive clusters
for clust = 1:length(pos_clust)
    
    %Find all significant voxels for current cluster across all times
    %[sig_grid_voxels, col] = find(clust_info.pos_clust_ids(:,:) == pos_clust(clust));
    
    %Find all significant voxels within a critical timeframe
    [sig_grid_voxels, col] = find(clust_info.pos_clust_ids(:,timeframe) == pos_clust(clust));

    %Aggergate voxels
    all_sig_pos_voxels = [all_sig_pos_voxels; sig_grid_voxels];
    
end

%Loop over postive clusters
for clust = 1:length(neg_clust)
    
    %Find all significant voxels for current cluster across all times
    %[sig_grid_voxels, col] = find(clust_info.neg_clust_ids(:,:) == neg_clust(clust));
    
    %Find all significant voxels within a critical timeframe
    [sig_grid_voxels, col] = find(clust_info.neg_clust_ids(:,timeframe) == neg_clust(clust));

    %Aggergate voxels
    all_sig_neg_voxels = [all_sig_neg_voxels; sig_grid_voxels];
    
end

%Find unique voxel positions
all_sig_pos_voxels = unique(all_sig_pos_voxels);
all_sig_neg_voxels = unique(all_sig_neg_voxels);

%Select positive clusters
if isequal(sig_clusters, 'p')
    
    %Directory name
    clust_dir_name = 'Positive Significant Voxels';
    
    %Combine pos neg cluster voxels
    all_sig_clust_voxels = all_sig_pos_voxels;

%Select negative clusters
elseif isequal(sig_clusters, 'n')
    
    %Directory name
    clust_dir_name = 'Negative Significant Voxels';
    
    %Combine pos neg cluster voxels
    all_sig_clust_voxels = all_sig_neg_voxels;

%Select positive and negative clusters 
elseif isequal(sig_clusters, 'b')
    
    %Directory name
    clust_dir_name = 'Positive and Negative Significant Voxels';
    
    %Combine pos neg cluster voxels
    all_sig_clust_voxels = unique([all_sig_pos_voxels; all_sig_neg_voxels]);

end

%Cut clustering data to only those voxels in significant neg and post clusters
clustering_matrix = clustering_matrix(all_sig_clust_voxels,:);

%If considering the grid voxels include their adjacent voxels
if isequal(voxel_method,1)
    
    %Select the significant voxel subset from matrix of adjacent voxels and
    %outInd position
    all_adjacent_voxels_cell = all_adjacent_voxels_cell(all_sig_clust_voxels);
    
end

%Only select significant voxels from grid_out
grid_outInd = grid_outInd(all_sig_clust_voxels);

%Rename for saving
kcluster_grid_outInd = grid_outInd;

%% Subcortical Atlas Info

%Subcortical ROI

%Load subcortical ROI voxel info - includes cerebellum
load(fullfile(subcortical_ROI_dir,'voxelInfo.mat'))

%Specify image size (should be 91x109x91);
sizeMat = [91,109,91]; 

%Get the list of index of useful voxels sizeMat(2:4), is the range of X, Y, Z indices. sizeMat(1) is time
subcort_outInd = mat2ind(voxelInfo,sizeMat); 

%Cortical ROI

%Load subcortical ROI voxel info - includes cerebellum
load(fullfile(cortical_ROI_dir,'voxelInfo.mat'))

%Specify image size (should be 91x109x91);
sizeMat = [91,109,91]; 

%Get the list of index of useful voxels sizeMat(2:4), is the range of X, Y, Z indices. sizeMat(1) is time
cort_outInd = mat2ind(voxelInfo,sizeMat); 

%% K-means Clustering and Plotting

%Loop over k-values
for cluster_num = 3%2:max_clusters %loops through k = 2 to the value of max_clusters
    
    %Define clusting save directory
    output_dir = fullfile(save_dir,ROI_name,clust_dir_name,['Cluster num ' num2str(cluster_num)]);
    
    %Make cluster directory if missing
    if ~exist(output_dir,'dir')
        
        mkdir(output_dir); 
        
    end
    
    %Run K-means
    [clusters, C, sumD, D] = kmeans(clustering_matrix(:, timeframe) , cluster_num, 'distance', 'correlation', ...
            'start', 'cluster', 'maxiter', maxiter, 'display', 'final', 'replicates', replicate_num);

    %Save k-means cluster info
    cd(output_dir)
    save(fullfile(output_dir,'clusterinfo.mat'), 'clusters', 'C', 'sumD', 'D', 'kcluster_grid_outInd');

    %% Plot Cluster Timecourse Whole Kmeans Cluster
    
    %Color
    color = {'r',[1 .54 0],'b',[0 .6 0],[1 .8 0],'c',[0.5 0 1],'m','k',[.5 .5 .5],[.93 .5 .73],...
        [0 0.5 1],[0.9 0.6 0.2], [0.9 0.1 0.6], [0.4 0.7 0.6], [0.2 0.1 0.5], [0.7 0.7 0.4], [0.1 0.5 0.3], 'g', [0.6 0.3 0.9], [0.7, 0.4, 0.4]};
    
    %Timepoint
    time = [-20:20];
    
    %Setup figure
    figure
    hold on
    
    title(['Cluster BOLD Timecourses K = ',num2str(cluster_num)])
    ylabel('Percent Change BOLD')
    xlabel('Time (s)')
    ylim([-0.15,  0.2])
    xlim([-20, 20])
    
    %Reference lines
    plot([0,0],[-0.15 ,0.2], 'k')
    plot([-20,20],[0,0], 'k')
    
    %Reset variables
    plot_handle = [];
    legendInfo = [];
    
    %Average BOLD over clusters
    for clust = 1:cluster_num
        
        %Plot timecourse of average BOLD signal within clusters
        plot_handle(clust) = plot(time,nanmean(clustering_matrix(clusters == clust,:),1),'Color', color{clust},'LineWidth',2);
        
        %Plot SEM (standard deviation across all voxels at each time point
        %divided by square root of the number of subjects (n = 99)
        plot(time,nanmean(clustering_matrix(clusters == clust,:),1) + std(clustering_matrix(clusters == clust,:),0,1)/sqrt(size(group_CP_BOLD_PC_data,3)),...
            'Color', color{clust},'LineWidth',1);
        plot(time,nanmean(clustering_matrix(clusters == clust,:),1) - std(clustering_matrix(clusters == clust,:),0,1)/sqrt(size(group_CP_BOLD_PC_data,3)),...
            'Color', color{clust},'LineWidth',1);
        
        %Store figure legend info
        legendInfo{clust} = ['Cluster ', num2str(clust), ' n = ', num2str(sum(clusters == clust))];
        
    end   
    
    %Plot legend
    lgd = legend(plot_handle,legendInfo, 'Location', 'Northwest');
    
    %Save figure
    cd(output_dir)
    savefig(['Cluster timecourses k = ', num2str(cluster_num)])
    
    close
    
    %% Plot Subcortical and Cortical Separated in K-means
    
    %Find subcortical voxels
    subcortical_voxels = ismember(grid_outInd,subcort_outInd);
    
    %Note: there are 15 voxels near subcortical sights that are not in
    %either cortical or subcortical ROIs - this line adds them to the
    %subcoritcal ROI just as has been implemented for the visualization of
    %the k-means voxels over time.
    subcortical_voxels = subcortical_voxels + (not(ismember(grid_outInd,subcort_outInd)) & not(ismember(grid_outInd,cort_outInd)));
    
    %Find cortical voxels
    cortical_voxels = ismember(grid_outInd,cort_outInd);
    
    %Color
    color = {'r',[1 .54 0],'b',[0 .6 0],[1 .8 0],'c',[0.5 0 1],'m','k',[.5 .5 .5],[.93 .5 .73],...
        [0 0.5 1],[0.9 0.6 0.2], [0.9 0.1 0.6], [0.4 0.7 0.6], [0.2 0.1 0.5], [0.7 0.7 0.4], [0.1 0.5 0.3], 'g', [0.6 0.3 0.9], [0.7, 0.4, 0.4]};
    
    %Timepoint
    time = [-20:20];
    
    %Setup figure
    figure
    hold on
    
    title(['Cluster BOLD Timecourses - Subcortical/Cortical - K = ',num2str(cluster_num)])
    ylabel('Percent Change BOLD')
    xlabel('Time (s)')
    ylim([-0.15,  0.2])
    xlim([-20, 20])
    
    %Reference lines
    plot([0,0],[-0.15 ,0.2], 'k')
    plot([-20,20],[0,0], 'k')
    
    %Reset variables
    plot_handle = [];
    legendInfo = [];

    %Average BOLD over clusters
    for clust = 1:cluster_num
             
        %Plot timecourse of average BOLD signal within clusters - Cortical
        %Voxels
        plot_handle(clust) = plot(time,nanmean(clustering_matrix(clusters == clust & cortical_voxels,:),1),'Color', color{clust},'LineWidth',2);
        
        %Plot SEM (standard deviation across all voxels at each time point
        %divided by square root of the number of subjects (n = 99)
        plot(time,nanmean(clustering_matrix(clusters == clust & cortical_voxels,:),1) + std(clustering_matrix(clusters == clust & cortical_voxels,:),0,1)/sqrt(size(group_CP_BOLD_PC_data,3)),...
            'Color', color{clust},'LineWidth',1);
        plot(time,nanmean(clustering_matrix(clusters == clust & cortical_voxels,:),1) - std(clustering_matrix(clusters == clust & cortical_voxels,:),0,1)/sqrt(size(group_CP_BOLD_PC_data,3)),...
            'Color', color{clust},'LineWidth',1);
        
        %Plot timecourse of average BOLD signal within clusters -
        %Subcortical Voxels
        plot_handle(clust) = plot(time,nanmean(clustering_matrix(clusters == clust & subcortical_voxels,:),1),'Color', color{clust},'LineWidth',2,'LineStyle','--');
        
        %Plot SEM (standard deviation across all voxels at each time point
        %divided by square root of the number of subjects (n = 99)
        plot(time,nanmean(clustering_matrix(clusters == clust & subcortical_voxels,:),1) + std(clustering_matrix(clusters == clust & subcortical_voxels,:),0,1)/sqrt(size(group_CP_BOLD_PC_data,3)),...
            'Color', color{clust},'LineWidth',1);
        plot(time,nanmean(clustering_matrix(clusters == clust & subcortical_voxels,:),1) - std(clustering_matrix(clusters == clust & subcortical_voxels,:),0,1)/sqrt(size(group_CP_BOLD_PC_data,3)),...
            'Color', color{clust},'LineWidth',1);
        
        %Store figure legend info
        legendInfo{clust} = ['Cluster ', num2str(clust), ' cortical n = ', num2str(sum(clusters == clust & cortical_voxels)),...
            ' subcortical n = ', num2str(sum(clusters == clust & subcortical_voxels))];
        
    end   
    
    %Plot legend
    lgd = legend(plot_handle,legendInfo, 'Location', 'Northwest');
    
    %Save figure
    cd(output_dir)
    savefig(['Cluster timecourses cortical subcortical k = ', num2str(cluster_num)])
    
    close
    
    %% Generate Nii and Mat Spatial Extent Map - All Voxels in Kmeans Cluster

    % load cluster .mat
    load(fullfile(output_dir,'clusterinfo.mat'));
    
    %Define new directory for cluster ROIs
    cluster_roi_dir = fullfile(output_dir, 'Cluster_roi');
    
    %Create cluster ROI if not made
    if ~exist(cluster_roi_dir,'dir')
        
        mkdir(cluster_roi_dir);
        
    end
    
    %Rename cluster label variable
    cluster_labels = clusters;
    
    %Find the total number of clusters
    num_clusters = max(cluster_labels);
    
    %Loop over clusters
    for clust = 1:num_clusters % want to create a .nii for each cluster
        
        %Find voxels in particular cluster
        curr_cluster = (cluster_labels == clust);
        
        %Total number of grid voxels array
        total_voxels = 1:size(curr_cluster,1);
        
        %If considering the grid voxels include their adjacent voxels 
        if isequal(voxel_method, 1)
        
            %Find adjacent voxels for curr_cluster voxels

            %Initialize array variable used to covert adjacent cell to array
            adj_voxels_array = [];

            %Loop over cluster voxels and find adjacent voxels
            for voxel = total_voxels(curr_cluster)

                %Store all adjacent voxels for each sig_grid_voxels in array
                adj_voxels_array = [adj_voxels_array; cell2mat(all_adjacent_voxels_cell(voxel))];

            end

            %METHOD 1 - GRID + ADJACENT VOXELS
            %Find xzy files in whole brain from subset of grid voxels that are
            %in a specific cluster; NOTE: the main matrix being cut down
            %via the curr cluster and grid/adjacent voxels must be the
            %total whole brain grey matter mask
            %curr_xyz = vXYZ(:,(ismember(outInd, [grid_outInd(curr_cluster); adj_voxels_array])));
            curr_xyz = temp_vXYZ(:,(ismember(outInd, [grid_outInd(curr_cluster); adj_voxels_array])));

        else
            
            %METHOD 2 - ALL VOXELS OR JUST GRID VOXELS 
            %Find xzy files in whole brain from subset of grid voxels that are
            %in a specific cluster
            %curr_xyz = vXYZ(:,(ismember(outInd, grid_outInd(curr_cluster))));
            curr_xyz = temp_vXYZ(:,(ismember(outInd, grid_outInd(curr_cluster))));

        end
        
        %Find grid voxels among total brain grey matter masked voxels
        %grid_to_whole_idx = ismember(outInd, grid_outInd);
    
        %maroi_pointlist is the Marsbar function to create an ROI
        act_roi = maroi_pointlist(struct('XYZ', curr_xyz, 'mat', V), 'vox'); 

        %Save ROI to MarsBaR ROI file
        cd(cluster_roi_dir)
        saveroi(act_roi, fullfile(cluster_roi_dir, ['cluster ' num2str(clust) ' of ' num2str(num_clusters) ' roi' '.mat']));
        save_as_image(act_roi, fullfile(cluster_roi_dir, ['cluster ' num2str(clust) ' of ' num2str(num_clusters) ' roi' '.nii']));

    end
    
    %% Generate Nii and Mat Spatial Extent Map - Cortical and Subcortical in Kmeans Cluster

    % load cluster .mat
    load(fullfile(output_dir,'clusterinfo.mat'));
    
    %Define new directory for cluster ROIs
    cluster_roi_dir = fullfile(output_dir, 'Cluster_roi');
    
    %Create cluster ROI if not made
    if ~exist(cluster_roi_dir,'dir')
        
        mkdir(cluster_roi_dir);
        
    end
    
    %Rename cluster label variable
    cluster_labels = clusters;
    
    %Find the total number of clusters
    num_clusters = max(cluster_labels);
    
    %Loop over cortical and subcortical voxels
    voxel_region = {'cortical','subcortical'};
    
    %Loop over regions
    for region = 1:length(voxel_region)
        
        %Select current region
        current_region = voxel_region{region};
    
        %Loop over clusters
        for clust = 1:num_clusters % want to create a .nii for each cluster

            %Find voxels in particular cluster
            curr_cluster = (cluster_labels == clust & eval([current_region,'_voxels']));

            %Total number of grid voxels array
            total_voxels = 1:size(curr_cluster,1);

            %If considering the grid voxels include their adjacent voxels 
            if isequal(voxel_method, 1)

                %Find adjacent voxels for curr_cluster voxels

                %Initialize array variable used to covert adjacent cell to array
                adj_voxels_array = [];

                %Loop over cluster voxels and find adjacent voxels
                for voxel = total_voxels(curr_cluster)

                    %Store all adjacent voxels for each sig_grid_voxels in array
                    adj_voxels_array = [adj_voxels_array; cell2mat(all_adjacent_voxels_cell(voxel))];

                end

                %METHOD 1 - GRID + ADJACENT VOXELS
                %Find xzy files in whole brain from subset of grid voxels that are
                %in a specific cluster; NOTE: the main matrix being cut down
                %via the curr cluster and grid/adjacent voxels must be the
                %total whole brain grey matter mask
                %curr_xyz = vXYZ(:,(ismember(outInd, [grid_outInd(curr_cluster); adj_voxels_array])));
                curr_xyz = temp_vXYZ(:,(ismember(outInd, [grid_outInd(curr_cluster); adj_voxels_array])));

            else

                %METHOD 2 - ALL VOXELS OR JUST GRID VOXELS 
                %Find xzy files in whole brain from subset of grid voxels that are
                %in a specific cluster
                %curr_xyz = vXYZ(:,(ismember(outInd, grid_outInd(curr_cluster))));
                curr_xyz = temp_vXYZ(:,(ismember(outInd, grid_outInd(curr_cluster))));

            end

            %Find grid voxels among total brain grey matter masked voxels
            %grid_to_whole_idx = ismember(outInd, grid_outInd);

            %maroi_pointlist is the Marsbar function to create an ROI
            act_roi = maroi_pointlist(struct('XYZ', curr_xyz, 'mat', V), 'vox'); 

            %Save ROI to MarsBaR ROI file
            cd(cluster_roi_dir)
            saveroi(act_roi, fullfile(cluster_roi_dir, ['cluster ' num2str(clust) ' of ' num2str(num_clusters) ' ' current_region ' roi' '.mat']));
            save_as_image(act_roi, fullfile(cluster_roi_dir, ['cluster ' num2str(clust) ' of ' num2str(num_clusters) ' ' current_region ' roi' '.nii']));

        end
        
    end
    
    %% Plot Silhouette Values
 
    %Color
    color = {'r',[1 .54 0],'b',[0 .6 0],[1 .8 0],'c',[0.5 0 1],'m','k',[.5 .5 .5],[.93 .5 .73],[0 0.5 1],[0.9 0.6 0.2], [0.9 0.1 0.6], [0.4 0.7 0.6], [0.2 0.1 0.5], [0.7 0.7 0.4], [0.1 0.5 0.3], 'g', [0.6 0.3 0.9], [0.7, 0.4, 0.4]};

    %Find silhouette value using the correlation distance (same distance
    %measure used for k-means clustering of these data).
    silh = silhouette(clustering_matrix(:,timeframe),clusters, 'correlation');

    %Plot
    figure
    hold on
    
    title(['K-means cluster Silhouette Values (k = ', num2str(cluster_num),')'])
    ylabel('Silouhette Value - Correlation')
    xlabel('Voxel - sorted by cluster and silouhette value')
    
    ylim([-0.2 0.8])

    %Setup bar_count variable
    bar_count_end = 0;

    %Average BOLD over clusters
    for clust = 1:cluster_num

        %Add 1 to bar count 
        bar_count_begin = bar_count_end + 1;
        
        %Define color
        if isstring(color{clust})
            
            current_color = color{clust};
            
        else
            
            current_color = cell2mat(color(clust));
            
        end

        %Plot timecourse of average BOLD signal within clusters
        plot_handle(clust) = bar([bar_count_begin:bar_count_end+sum(clusters == clust)],sort(silh(clusters == clust,:),'descend'),1,'FaceColor',current_color);

        %Plot the mean silhouette value marker
        plot([(bar_count_begin + bar_count_end + sum(clusters == clust))/2], mean(silh(clusters == clust,:)),'ok')

        %Add total number of voxels for a cluster to bar_count
        bar_count_end = bar_count_end + sum(clusters == clust);

        %Store figure legend info
        legendInfo{clust} = ['Cluster ', num2str(clust), ' n = ', num2str(sum(clusters == clust))];

    end   

    %Save figure
    cd(output_dir)
    savefig(['K-means clusters silhouette figure.fig'])
    
    close
    
end

%% Obtain silhouette values

silhouette_dir = fullfile(output_base, 'silhouette_values');

if ~exist(silhouette_dir,'dir')
    mkdir(silhouette_dir);
end

k_vals = 2:max_clusters;

all_avg_cluster_avg_values = zeros(length(k_vals),1);
all_avg_avg_values = NaN(length(k_vals),max(k_vals));

shift = k_vals(1) - 1;
for curr_cluster = k_vals
    disp(['Curr cluster is ',num2str(curr_cluster)])
    
    % load the data
    load(fullfile(output_base,['Cluster_num ' num2str(curr_cluster)],'clusterinfo.mat'));
        
    avg_val_per_cluster = NaN(curr_cluster,1);
    
    for j = 1:curr_cluster
        indices = find(clusters==j);
        
        % compute the silhouette value
        [s, h] = silhouette(final_data, clusters, 'correlation');
        
        %take average sil value within cluster
        avg_val_per_cluster(j) = mean(s(indices));
        all_avg_avg_values(curr_cluster-shift,j) = avg_val_per_cluster(j);
    end
    all_avg_cluster_avg_values(curr_cluster-shift) = mean(avg_val_per_cluster);
    save(fullfile(silhouette_dir,'silhouette.mat'),'all_avg_avg_values');
end

%% Plotting silhouette values

% x_range = size(all_avg_avg_values,1); % plot all values
x_range = MAX_clusters; % max value

% plotting all values
h = figure;
hold on
for a = 1:x_range-1
    scatter(a+1*ones(1,a+1),all_avg_avg_values(a, 1:a+1), 'filled')
end
xlabel('Number of clusters')
ylabel('Average silhouette value per cluster')
title('All cluster average silhouette values')
axis([1.5 MAX_clusters+0.5 0 0.7])
% saveas(h, fullfile(silhouette_dir,['all_values_2_',num2str(x_range),'.tif']))
% close(h)

% plotting all values, connecting lowest value
h = figure;
hold on
for a = 1:x_range-1
    scatter(a+1*ones(1,a+1),all_avg_avg_values(a, 1:a+1), 'filled')
    lowest_val(a) = min(all_avg_avg_values(a,1:a+1));
end
plot(2:x_range, lowest_val, 'k')
xlabel('Number of clusters')
ylabel('Average silhouette value per cluster')
title('All cluster average silhouette values, lowest value connected')
axis([1.5 MAX_clusters+0.5 0 0.7])
% saveas(h, fullfile(silhouette_dir,['all_values_lowest_connected_2_',num2str(x_range),'.tif']))
% close(h)

% plotting only lowest value
h = figure;
hold on
for a = 1:x_range-1
    lowest_val(a) = min(all_avg_avg_values(a,1:a+1));
    scatter(a+1,lowest_val(a), 'filled','b')
end
xlabel('Number of clusters')
ylabel('Minimum average silhouette value per cluster')
title('Minimum cluster average silhouette values')
axis([1.5 MAX_clusters+0.5 0 0.55])
% saveas(h, fullfile(silhouette_dir,['lowest_values_2_',num2str(x_range),'.tif']))
% close(h)

% save variables
% save(fullfile(silhouette_dir,'cluster_silhouette_averages.mat'),'all_avg_avg_values','all_avg_cluster_avg_values')
