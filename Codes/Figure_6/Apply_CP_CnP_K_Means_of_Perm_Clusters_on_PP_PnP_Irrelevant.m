%% Apply K-means clusters from CP and CnP data to PP and PnP data

%This code will take the k-means clusters from the reported perceived and
%not perceived trial class (CP and CnP) and apply them onto the PP and PnP
%data. 

%Note: data needs to be cut down to the CP and CnP statisitcally sig voxels
%because these are the voxels used in k-means clustering

% Name: Sharif Kronemer
% Date: 6/5/2021
% Modified: 6/24/2021

clear

%% Prompts

%Run location
% prompt1 = 'Run location (s/l): ';
% run_location = input(prompt1, 's');
run_location = 's';

%Combined Kmeans clusters and NRP kmeans clusters
prompt1 = 'Select kmeans clusters [NRP_CP_CnP, Combined_CP_CnP]: ';
source_clusters = input(prompt1, 's');

%Define name of source cluster
if isequal(source_clusters,'NRP_CP_CnP')
    
    source_filename = 'NRP Task CP CnP Kclust';
    
elseif isequal(source_clusters,'Combined_CP_CnP')
    
    source_filename = 'Combined Task CP CnP Kclust';
    
end

%Select Positive or Negative Clusters
prompt2 = 'Positive Clusters, Negative Clusters, or Both (p,n,b): ';
sig_clusters = input(prompt2, 's');

%Stat sig irrelevant voxels only
prompt3 = 'Sig irrelevant voxels only (y,n): ';
sig_voxels_only = input(prompt3, 's');

%% Parameters 

%Variable used for establishing directories and folder names
ROI_name = 'Whole Brain';

%Condience Score
confidence_score = 0.75;

%Number of permutations
num_permutation = 5000;

%Center and quadrant name variable 
location_name = 'cent_quad';

%% Add Paths and Directories

%Save folder
if isequal(sig_voxels_only,'n')
    
    save_folder = ['PP_PnP_',location_name,'_score_',num2str(confidence_score),'_',num2str(num_permutation),'_allvoxels'];

elseif isequal(sig_voxels_only,'y')
    
    save_folder = ['PP_PnP_',location_name,'_score_',num2str(confidence_score),'_',num2str(num_permutation)];
    
end
    
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
    template_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Subject Analysis/191AM/MRI Analysis/Extracted voxel data/Run 1 Movie/axialSlices'; 

    %Define CP/CnP Stats and k clusters
    if isequal(source_clusters,'NRP_CP_CnP')

        %CP CnP perm stat cluster directory
        CP_CnP_stat_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/MRI Analysis/Percent Change Permutation Stats',...
            [ROI_name,'_CP_minus_CnP_relevant_',location_name,'_threshold_perm_',num2str(num_permutation)]);   

        %CP/CnP clusters
        CP_CnP_kclusters = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/MRI Analysis/K means Clusters/CP_CnP_trials/Whole Brain';

    elseif isequal(source_clusters,'Combined_CP_CnP')

        %CP CnP perm stat cluster directory
        CP_CnP_stat_dir = fullfile('/mnt/Data8/HNCT NRP and Report Only Paradigm Analyses/Analysis/MRI Analysis/Percent Change Permutation Stats/Whole Brain CP minus CnP Avg 5000');
        
        %CP/CnP clusters
        CP_CnP_kclusters = '/mnt/Data8/HNCT NRP and Report Only Paradigm Analyses/Analysis/MRI Analysis/K means Clusters/Whole Brain';
    
    end
    
    %PP PnP perm stat cluster directory
    PP_PnP_stat_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/MRI Analysis/Percent Change Permutation Stats',...
        [ROI_name,'_PP_minus_PnP_irrelevant_',location_name,'_threshold_score_thres_',num2str(confidence_score),'_perm_',num2str(num_permutation)]);
   
    %MNI brain
    ROI_dir = ['/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Permutation Statistical Analysis/Brain Voxel Subset/', ROI_name]; %Cortical/Subcortical voxels
    
    %PP and PnP Data dir
    data_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Data/MRI/BOLD PC';
    
    %Save dir
    save_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/MRI Analysis/K means Clusters/PP_PnP_trials/Whole Brain',source_filename,save_folder);
    
    %Cortical and Subcortical ROI (Custom made cortical and subcortical ROI from MarsBaR and harvard atlas)
    subcortical_ROI_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/Subcortical ROI';
    cortical_ROI_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/Marsbar ROIs/Cortical ROI';
    
end

%Make save folder
mkdir(save_dir)

%% Parameters

% %Clustering parameters
max_clusters = 10;  %max number of clusters, these many clusters will be attempted
timeframe = 21:31; %the epoch in the precent change map data that will be clustered (e.g., 21s to 36s represents the time from face onset to the end of the 15s post-stim window)

%% Generate Group Percent Change Variable [time x voxel x subjects]

disp('Load BOLD PC Data')

%Load PP and PnP PC Data
cd(data_dir)
load(['Group_PC_BOLD_irrelevant_PP_PnP_',location_name,'_threshold_score_thres_',num2str(confidence_score),'_data.mat'])

%Subtract main data from subtract data
avg_group_pc_data = nanmean(group_PP_BOLD_PC_data,3) - nanmean(group_PnP_BOLD_PC_data,3);
        
%% Setup Image Template

%Voxel Method defines if all or a subset of voxels will be considered in clustering
%0 = use all voxels 
%1 = use grid/ROI subset voxels
voxel_method = 1;

%Load ROI voxelInfo
cd(ROI_dir)
load(fullfile(ROI_dir,'voxelInfo.mat'));

%{
%Define XYZ position and select grid subset
try
    
    %Voxel info
    vXYZ = xyz.vXYZ;
    
    %Image matrix
    %V = xyz.mat;
    
catch
    
    %Voxel info
    vXYZ = voxelInfo;
    
    %Image matrix
    %V = transform_matrix;
    
end
%}

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

%% Load Permutation Clusters and Select Significant Clusters (p < 0.05) - Permutation clusters coming from CP-CnP testing

%Note: CP/CnP stats were used to determine the voxels for CP/CnP k-means
%clustering and therefore these stats need to be applied to the PP/PnP data

%Load permutation cluster information
cd(CP_CnP_stat_dir)
load('pc_cluster_sumt_rand_005.mat')

%Find significant clusters pvalue < 0.05 - Negative and positive clusters independently
pos_clust = find(clust_info.pos_clust_pval < 0.05);
neg_clust = find(clust_info.neg_clust_pval < 0.05);

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

% %Find unique voxel positions
% all_sig_pos_voxels = unique(all_sig_pos_voxels);
% all_sig_neg_voxels = unique(all_sig_neg_voxels);

%Select positive clusters
if isequal(sig_clusters, 'p')
    
    %Directory name
    clust_dir_name = 'Pos Sig Voxels';
    
    %Combine pos neg cluster voxels
    rel_all_sig_clust_voxels = unique(all_sig_pos_voxels);

%Select negative clusters
elseif isequal(sig_clusters, 'n')
    
    %Directory name
    clust_dir_name = 'Neg Sig Voxels';
    
    %Combine pos neg cluster voxels
    rel_all_sig_clust_voxels = unique(all_sig_neg_voxels);

%Select positive and negative clusters 
elseif isequal(sig_clusters, 'b')
    
    %Directory name
    clust_dir_name = 'Pos Neg Sig Voxels';
    clust_dir_full_name = 'Positive and Negative Significant Voxels';
    
    %Combine pos neg cluster voxels
    rel_all_sig_clust_voxels = unique([all_sig_pos_voxels; all_sig_neg_voxels]);

end

%Cut clustering data to only those voxels in significant neg and post clusters
clustering_matrix = clustering_matrix(rel_all_sig_clust_voxels,:);

%If considering the grid voxels include their adjacent voxels
if isequal(voxel_method,1)
    
    %Select the significant voxel subset from matrix of adjacent voxels and
    %outInd position
    all_adjacent_voxels_cell = all_adjacent_voxels_cell(rel_all_sig_clust_voxels);
    
end

%Only select significant voxels from grid_out
grid_outInd = grid_outInd(rel_all_sig_clust_voxels);

%% Load Permutation Clusters and Select Significant Clusters (p < 0.05) - Permutation clusters coming from PP-PnP testing

%Select the voxels to use for clustering
if isequal(sig_voxels_only,'y')

    %Load permutation cluster information
    cd(PP_PnP_stat_dir)
    load('pc_cluster_sumt_rand_005.mat')

    %Find significant clusters pvalue < 0.05 - Negative and positive clusters independently
    pos_clust = find(clust_info.pos_clust_pval < 0.05);
    neg_clust = find(clust_info.neg_clust_pval < 0.05);

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

    %Select positive clusters
    if isequal(sig_clusters, 'p')

        %Find unique voxel positions
        irrel_all_sig_clust_voxels = unique(all_sig_pos_voxels);

    %Select negative clusters
    elseif isequal(sig_clusters, 'n')

        %Find unique voxel positions
        irrel_all_sig_clust_voxels = unique(all_sig_neg_voxels);

    %Select positive and negative clusters 
    elseif isequal(sig_clusters, 'b')

        %Combine pos neg cluster voxels
        irrel_all_sig_clust_voxels = unique([all_sig_pos_voxels; all_sig_neg_voxels]);

    end

    %Find common sig voxels and use only the CP/CnP stat sig voxels that
    %are also sig for the PP/PnP contrast
    union_perm_stat_voxels = ismember(rel_all_sig_clust_voxels, irrel_all_sig_clust_voxels);
    
end

%% Cortical and Subcortical Atlas Info

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
for cluster_num = 1:max_clusters %loops through k = 2 to the value of max_clusters
    
    %Define clusting save directory
    output_dir = fullfile(save_dir,clust_dir_name,['Cluster num ' num2str(cluster_num)]);
    mkdir(output_dir)
    
    %Load k-means clustering info CP/CnP 
    cluster_info_dir = fullfile(CP_CnP_kclusters,clust_dir_full_name,['Cluster num ' num2str(cluster_num)]);
    %mkdir(cluster_info_dir)
    
    cd(cluster_info_dir)
    load('clusterinfo.mat')
    
    %For sig voxels only condition
    if isequal(sig_voxels_only,'y')
        
        %Check the number of voxels is equal between indices
        if not(isequal(size(clusters,1),size(union_perm_stat_voxels,1)))

           error('Indices do not match!')

        end

    end

    %% Plot Cluster Timecourse Whole Kmeans Cluster
    
    %Color
    color = {'r',[1 .54 0],'b',[0 .6 0],[1 .8 0],'c',[0.5 0 1],'m','k',[.5 .5 .5],[.93 .5 .73],...
        [0 0.5 1],[0.9 0.6 0.2], [0.9 0.1 0.6], [0.4 0.7 0.6], [0.2 0.1 0.5], [0.7 0.7 0.4],...
        [0.1 0.5 0.3], 'g', [0.6 0.3 0.9], [0.7, 0.4, 0.4]};
    
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
    for kclust = 1:cluster_num
        
        %Plot timecourse of average BOLD signal within clusters
        if isequal(sig_voxels_only,'n')
            
            %Consider all voxels in that k cluster
            plot_handle(kclust) = plot(time,nanmean(clustering_matrix(clusters == kclust,:),1),'Color', color{kclust},'LineWidth',2);
            
            %Plot SEM (standard deviation across all voxels at each time point
            %divided by square root of the number of subjects 
            plot(time,nanmean(clustering_matrix(clusters == kclust,:),1) + std(clustering_matrix(clusters == kclust,:),0,1)/sqrt(size(group_PP_BOLD_PC_data,3)),...
                'Color', color{kclust},'LineWidth',1);
            plot(time,nanmean(clustering_matrix(clusters == kclust,:),1) - std(clustering_matrix(clusters == kclust,:),0,1)/sqrt(size(group_PP_BOLD_PC_data,3)),...
                'Color', color{kclust},'LineWidth',1);

            %Legend
            legendInfo{kclust} = ['Cluster ', num2str(kclust), ' n = ', num2str(sum(clusters == kclust))];

        elseif isequal(sig_voxels_only,'y')
            
            %Consider only voxels that are statistically sig for PP/PnP in
            %that k cluster
            plot_handle(kclust) = plot(time,nanmean(clustering_matrix(clusters == kclust & union_perm_stat_voxels,:),1),'Color', color{kclust},'LineWidth',2); %Only kmeans voxels with common perm stat voxels
            
            %Plot SEM (standard deviation across all voxels at each time point
            %divided by square root of the number of subjects 
            plot(time,nanmean(clustering_matrix(clusters == kclust & union_perm_stat_voxels,:),1) + std(clustering_matrix(clusters == kclust & union_perm_stat_voxels,:),0,1)/sqrt(size(group_PP_BOLD_PC_data,3)),...
                'Color', color{kclust},'LineWidth',1);
            plot(time,nanmean(clustering_matrix(clusters == kclust & union_perm_stat_voxels,:),1) - std(clustering_matrix(clusters == kclust & union_perm_stat_voxels,:),0,1)/sqrt(size(group_PP_BOLD_PC_data,3)),...
                'Color', color{kclust},'LineWidth',1);
                      
            %Legend
            legendInfo{kclust} = ['Cluster ', num2str(kclust), ' n = ', num2str(sum(clusters == kclust & union_perm_stat_voxels))]; %Only count voxels with significance in common

        end
    
    end   
    
    %Plot legend
    lgd = legend(plot_handle,legendInfo, 'Location', 'Northwest');
    
    %Save figure
    cd(output_dir)
    savefig(['Cluster timecourses k = ', num2str(cluster_num)])
    
    close
    
    %% Plot Subcortical and Cortical Separated in K-means
    
    %Find subcortical voxels
    subcort_voxels = ismember(grid_outInd,subcort_outInd);
    
    %Note: there are 15 voxels near subcortical sights that are not in
    %either cortical or subcortical ROIs - this line adds them to the
    %subcoritcal ROI just as has been implemented for the visualization of
    %the k-means voxels over time.
    subcort_voxels = subcort_voxels + (not(ismember(grid_outInd,subcort_outInd)) & not(ismember(grid_outInd,cort_outInd)));
    
    %Find cortical voxels
    cort_voxels = ismember(grid_outInd,cort_outInd);
    
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
        
        %Plot timecourse of average BOLD signal within clusters
        if isequal(sig_voxels_only,'y')

            %Plot timecourse of average BOLD signal within clusters - Cortical
            %Voxels
            plot_handle(clust) = plot(time,nanmean(clustering_matrix(clusters == clust & cort_voxels & union_perm_stat_voxels,:),1),'Color', color{clust},'LineWidth',2);

            %Plot SEM (standard deviation across all voxels at each time point
            %divided by square root of the number of subjects (n = 99)
            plot(time,nanmean(clustering_matrix(clusters == clust & cort_voxels & union_perm_stat_voxels,:),1) + std(clustering_matrix(clusters == clust & cort_voxels & union_perm_stat_voxels,:),0,1)/sqrt(size(group_PP_BOLD_PC_data,3)),...
                'Color', color{clust},'LineWidth',1);
            plot(time,nanmean(clustering_matrix(clusters == clust & cort_voxels & union_perm_stat_voxels,:),1) - std(clustering_matrix(clusters == clust & cort_voxels & union_perm_stat_voxels,:),0,1)/sqrt(size(group_PP_BOLD_PC_data,3)),...
                'Color', color{clust},'LineWidth',1);

            %Plot timecourse of average BOLD signal within clusters -
            %Subcortical Voxels
            plot_handle(clust) = plot(time,nanmean(clustering_matrix(clusters == clust & subcort_voxels & union_perm_stat_voxels,:),1),'Color', color{clust},'LineWidth',2,'LineStyle','--');

            %Plot SEM (standard deviation across all voxels at each time point
            %divided by square root of the number of subjects (n = 99)
            plot(time,nanmean(clustering_matrix(clusters == clust & subcort_voxels & union_perm_stat_voxels,:),1) + std(clustering_matrix(clusters == clust & subcort_voxels & union_perm_stat_voxels,:),0,1)/sqrt(size(group_PP_BOLD_PC_data,3)),...
                'Color', color{clust},'LineWidth',1);
            plot(time,nanmean(clustering_matrix(clusters == clust & subcort_voxels & union_perm_stat_voxels,:),1) - std(clustering_matrix(clusters == clust & subcort_voxels & union_perm_stat_voxels,:),0,1)/sqrt(size(group_PP_BOLD_PC_data,3)),...
                'Color', color{clust},'LineWidth',1);

            %Store figure legend info
            legendInfo{clust} = ['Cluster ', num2str(clust), ' cortical n = ', num2str(sum(clusters == clust & cort_voxels & union_perm_stat_voxels)),...
                ' subcortical n = ', num2str(sum(clusters == clust & subcort_voxels & union_perm_stat_voxels))];
        
        elseif isequal(sig_voxels_only,'n')
           
            %Plot timecourse of average BOLD signal within clusters - Cortical
            %Voxels
            plot_handle(clust) = plot(time,nanmean(clustering_matrix(clusters == clust & cort_voxels,:),1),'Color', color{clust},'LineWidth',2);

            %Plot SEM (standard deviation across all voxels at each time point
            %divided by square root of the number of subjects (n = 99)
            plot(time,nanmean(clustering_matrix(clusters == clust & cort_voxels,:),1) + std(clustering_matrix(clusters == clust & cort_voxels,:),0,1)/sqrt(size(group_PP_BOLD_PC_data,3)),...
                'Color', color{clust},'LineWidth',1);
            plot(time,nanmean(clustering_matrix(clusters == clust & cort_voxels,:),1) - std(clustering_matrix(clusters == clust & cort_voxels,:),0,1)/sqrt(size(group_PP_BOLD_PC_data,3)),...
                'Color', color{clust},'LineWidth',1);

            %Plot timecourse of average BOLD signal within clusters -
            %Subcortical Voxels
            plot_handle(clust) = plot(time,nanmean(clustering_matrix(clusters == clust & subcort_voxels,:),1),'Color', color{clust},'LineWidth',2,'LineStyle','--');

            %Plot SEM (standard deviation across all voxels at each time point
            %divided by square root of the number of subjects (n = 99)
            plot(time,nanmean(clustering_matrix(clusters == clust & subcort_voxels,:),1) + std(clustering_matrix(clusters == clust & subcort_voxels,:),0,1)/sqrt(size(group_PP_BOLD_PC_data,3)),...
                'Color', color{clust},'LineWidth',1);
            plot(time,nanmean(clustering_matrix(clusters == clust & subcort_voxels,:),1) - std(clustering_matrix(clusters == clust & subcort_voxels,:),0,1)/sqrt(size(group_PP_BOLD_PC_data,3)),...
                'Color', color{clust},'LineWidth',1);

            %Store figure legend info
            legendInfo{clust} = ['Cluster ', num2str(clust), ' cortical n = ', num2str(sum(clusters == clust & cort_voxels)),...
                ' subcortical n = ', num2str(sum(clusters == clust & subcort_voxels))];
                                  
        end
        
    end   
    
    %Plot legend
    lgd = legend(plot_handle,legendInfo, 'Location', 'Northwest');
    
    %Save figure
    cd(output_dir)
    savefig(['Cluster timecourses cortical subcortical k = ', num2str(cluster_num)])
    
    close
    
    %% Generate Nii and Mat Spatial Extent Map - All Voxels

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
    for kclust = 1:num_clusters % want to create a .nii for each cluster
        
        %Get the current cluster voxels
        if isequal(sig_voxels_only,'y')
            
            %Find voxels in particular cluster - only sig voxels
            curr_cluster = (cluster_labels == kclust & union_perm_stat_voxels);        
        
        else
            
            %Find voxels in particular cluster - all voxels
            curr_cluster = (cluster_labels == kclust);
        
        end
        
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
        
        %Save name different for sig vs all voxels
        if isequal(sig_voxels_only,'y')
            
            saveroi(act_roi, fullfile(cluster_roi_dir, ['cluster ' num2str(kclust) ' of ' num2str(num_clusters) ' sigv.mat']));
            save_as_image(act_roi, fullfile(cluster_roi_dir, ['cluster ' num2str(kclust) ' of ' num2str(num_clusters) ' sigv.nii']));

        else
            
          saveroi(act_roi, fullfile(cluster_roi_dir, ['cluster ' num2str(kclust) ' of ' num2str(num_clusters) '.mat']));
          save_as_image(act_roi, fullfile(cluster_roi_dir, ['cluster ' num2str(kclust) ' of ' num2str(num_clusters) '.nii']));

        end
        
    end
   
    %% Generate Nii and Mat Spatial Extent Map - Cortical and Subcortical in kmeans Cluster

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
    voxel_region = {'cort','subcort'};
    
    %Loop over regions
    for region = 1:length(voxel_region)
        
        %Select current region
        current_region = voxel_region{region};
    
        %Loop over clusters
        for kclust = 1:num_clusters % want to create a .nii for each cluster

            %Find voxels in particular cluster
            curr_cluster = (cluster_labels == kclust & eval([current_region,'_voxels']));

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
            saveroi(act_roi, fullfile(cluster_roi_dir, ['cluster ' num2str(kclust) ' of ' num2str(num_clusters) ' ' current_region '.mat']));
            save_as_image(act_roi, fullfile(cluster_roi_dir, ['cluster ' num2str(kclust) ' of ' num2str(num_clusters) ' ' current_region '.nii']));

        end
    
    end
    
    %% Plot Silhouette Values
 
    %Color
    color = {'r',[1 .54 0],'b',[0 .6 0],[1 .8 0],'c',[0.5 0 1],'m','k',[.5 .5 .5],[.93 .5 .73],[0 0.5 1],[0.9 0.6 0.2], ...
        [0.9 0.1 0.6], [0.4 0.7 0.6], [0.2 0.1 0.5], [0.7 0.7 0.4], [0.1 0.5 0.3], 'g', [0.6 0.3 0.9], [0.7, 0.4, 0.4]};

    %Define silhouette clusters
    if isequal(sig_voxels_only,'n')
        
        %Select all voxels
        silh_clusters = clusters;
        
    elseif isequal(sig_voxels_only,'y')
        
        %Select only the voxels that are sig for PP/PnP
        silh_clusters = clusters(union_perm_stat_voxels);

    end
    
    %Find silhouette value using the correlation distance (same distance
    %measure used for k-means clustering of these data).
    silh = silhouette(clustering_matrix(union_perm_stat_voxels,timeframe),silh_clusters, 'correlation');
    silh = silhouette(clustering_matrix(union_perm_stat_voxels,27:31),silh_clusters, 'correlation');

    %Plot silhouette
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

        %Plot bar of silh values
        plot_handle(clust) = bar([bar_count_begin:bar_count_end+sum(silh_clusters == clust)],sort(silh(silh_clusters == clust,:),'descend'),1,'FaceColor',current_color);
        %plot_handle(clust) = plot([bar_count_begin:bar_count_end+sum(silh_clusters == clust)],sort(silh(silh_clusters == clust,:),'descend'),current_color);

        %Plot the mean silhouette value marker
        plot([(bar_count_begin + bar_count_end + sum(silh_clusters == clust))/2], mean(silh(silh_clusters == clust,:)),'ok')

        %Add total number of voxels for a cluster to bar_count
        bar_count_end = bar_count_end + sum(silh_clusters == clust);

        %Store figure legend info
        legendInfo{clust} = ['Cluster ', num2str(clust), ' n = ', num2str(sum(silh_clusters == clust))];

    end   

    %Save figure
    cd(output_dir)
    savefig(['K-means clusters silhouette figure.fig'])
    
    close
    
end
