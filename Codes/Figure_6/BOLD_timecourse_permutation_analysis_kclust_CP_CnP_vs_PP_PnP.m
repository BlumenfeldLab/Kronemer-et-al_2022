%% MR Voxel Based Permutation Analysis - ROI Timecourse - CP/CnP vs PP/PnP

%This code variant will cut each ROI selected for timecourse analysis and
%plotting into the DAS, TPN, and DMN k-means cluster networks. 

%The purpose of this code is to take voxels from an ROI, average over those
%voxels (thus collapsing the space dimension), and find statically
%signficiant changes in time via the cluster based permutation method. This
%is completed both by CP vs CnP and CP-CnP vs baseline (prestimulus period)

%Note about ROI entry: experimenter must type in the ROI and atlas name in
%the code prior to running. 

%This code will:
%(1) Load percent change data and aggregate across subjects
%(2) Cut voxels to ROI size
%(3) Cluster based permutation analysis 
%(4) Generate plots ROI timecourse plots with significant time periods
%noted

%Written by: Sharif I. Kronemer
%Date: 9/3/2020
%Modified: 7/22/2021

clear

%% Prompts

% %Run location
% prompt1 = 'Run location (s/l): ';
% run_location = input(prompt1, 's');
run_location = 's'; %Hard code run location to server

%Select atlas
prompt2 = 'Atlas  [MarsBaR, Stanford, Xiao, Neudorfer, Custom]: ';
Atlas_name = input(prompt2, 's');

%Select ROIs
if isequal(Atlas_name, 'Stanford')
    
    %ROI cell
    ROI_cell= {'medial_parietal_Precuneus_02','motor_cortex_left_Sensorimotor_01','precuneus_medial_PFC_dorsal_DMN_01_04','operculum_FEF_Visuospatial_03_07',...
        'SMA_AC_anterior_Salience_03','inferior_parietal_gyrus_LECN_03_RECN_03','frontal_gyrus_LECN_01_RECN_01','medial_PFC_dorsal_DMN_01','precuneus_dorsal_DMN_04','middle_frontal_gyrus_LECN_02_RECN_02'};
    
        %{'frontal_cortex_operculum_LECN_01_RECN_01_Visuospatial_03_07','inferior_parietal_gyrus_Visuospatial_02_06', 'inferior_parietal_gyrus_ECN_Visuospatial','thalamus_caudate_putamen_Basal_Ganglia_01_02',...
        %'cerebellum_LECN_05_RECN_05',...
        %'middle_frontal_gyrus_LECN_02_RECN_02','V1_prim_Visual_01',...
        %'insula_anterior_Salience_02_05','anterior_Salience','Basal_Ganglia','dorsal_DMN','ventral_DMN','Sensorimotor','high_Visual',...
        %'post_Salience','prim_Visual','Visuospatial','LECN','RECN','right_left_ECN','prim_high_Visual','ant_post_Salience'};

elseif isequal(Atlas_name,'MarsBaR')
    
    %ROI cell
    ROI_cell = {'MNI_Frontal_Sup_Orb_R+L_roi', 'MNI_Fusiform_L+R_roi', 'MNI_Parietal_Sup_L+R_roi','MNI_Insula_L+R_roi','Thalamus L+R',...
        'Midbrain_Reticular_Formation_roi','MNI_Vermis_Whole_roi','MNI_Cerebelum_Crus1_Crus2_L+R_roi','MNI_Temporal_Pole_Mid_Med_L+R_roi'};
    
    %{'MNI_Parietal_Sup_L+R_roi'};%{'MNI_Temporal_Mid_R+L_roi','MNI_Hippocampus_ParaHippo_R+L_roi','MNI_Lateral_Occipital_L+R_roi','MNI_Temporal_L+R_roi_roi',};%{'MNI_Cingulate_L+R',...
%        'Midbrain_Reticular_Formation_roi','Intralaminar_Thalamic_Nuclei_v5_symetric_roi','MNI_Insula_L+R_roi','Thalamus L+R',...
%        'SMA L+R', 'DMNv2 Region 2','MNI_Cerebelum_Crus1_L+R_roi','MNI_Cerebelum_Crus2_L+R_roi','MNI_Frontal_Mid_L+R_roi',...
%        'MNI_Parietal_Inf_Sup_L+R','MNI_Striatum_L+R_roi','MNI_Vermis_Whole_roi'};
%   {'MNI_Amygdala_L+R_roi','MNI_Frontal_Cortex_L+R','MNI_Calcarine_L+R_roi','MNI_Fusiform_L+R_roi','MNI_Cingulate_L+R',...
%         'Midbrain_Reticular_Formation_roi','Intralaminar_Thalamic_Nuclei_v5_symetric_roi','MNI_Insula_L+R_roi','Thalamus L+R',...
%         'SMA L+R', 'DMNv2 Region 2','MNI_Cerebelum_Crus1_L+R_roi','MNI_Cerebelum_Crus2_L+R_roi','MNI_Frontal_Mid_L+R_roi',...
%         'MNI_Parietal_Inf_Sup_L+R','MNI_Striatum_L+R_roi','MNI_Vermis_Whole_roi'};

elseif isequal(Atlas_name,'Xiao')
    
    %ROI cell
    ROI_cell = {'Claustrum','Subthalamus'};
    
elseif isequal(Atlas_name,'Neudorfer')
    
    %ROI cell
    ROI_cell = {'Nucleus Basalis'};
    
elseif isequal(Atlas_name,'Custom')
    
    %ROI cell
    ROI_cell = {'NA_custom_ROI_L+R_roi'};

    
end

%% Parameters and Variables Names

%Location set
location_data_name = 'cent_quad';
location_folder_name = 'Center and Quadrant Irrelevant';

%Confidence score
confidence_score = 0.75;

%Are samples dependent (default is true)
dependent_samples = 'true';

%Define alpha threshold
p_threshold = 0.05;

%Define the number of permutations
num_permutations = 5000;

%Two-sided
two_sided = 'true';

%% Add Functions to Path

if isequal(run_location, 'l')
    
elseif isequal(run_location, 's')
    
    %Add SPM and supplementary functions to path
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Permutation Statistical Analysis/Permutation functions')
    addpath('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Time courses/Supplementary functions/fMRI Analysis/Code for Percent Change Maps/subRoutines')
    addpath('/mnt/Data5/Sharif Data F/HNCT fMRI Study/spm12_radiological')
    addpath(genpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Extract voxel data/Supplementary functions'));
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Permutation Statistical Analysis')

    %Save directory 
    save_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/MRI Analysis/Percent Change Timecourses/Permutation Statistical Analysis',...
        location_folder_name,['CP_CnP_vs_PP_PnP_trials_kcluster_confidence_score_',num2str(confidence_score)],[Atlas_name,'_ROIs']);
    
    %Combined data direcotry
    irrel_group_data_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Data/MRI/BOLD PC';
    rel_group_data_dir = '/mnt/Data8/HNCT NRP and Report Only Paradigm Analyses/Group Data/MRI/BOLD PC';
    
end

%% Loop over ROI cell

for roi = 1:length(ROI_cell)
       
    %Select current ROI
    ROI_name = ROI_cell{roi};
    
    disp(['Running ROI ', ROI_name])
        
    %ROI folder name
    ROI_folder = [ROI_name, ' Timecourse Perm ',num2str(num_permutations)];

    %Figure folder
    fig_folder = fullfile(save_dir, ROI_folder);
    mkdir(fig_folder)

    %MarsBaR Atlas
    if isequal(Atlas_name, 'MarsBaR')
        
        ROI_dir = fullfile('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/Marsbar ROIs',ROI_name); 

   %Stanford Atlas
   elseif isequal(Atlas_name, 'Stanford')
       
        %Atlas directory
        ROI_dir = fullfile('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/Stanford Functional ROIs',ROI_name);       
    
   %Xiao Atlas
   elseif isequal(Atlas_name, 'Xiao')
        
        %Atlas directory
        ROI_dir = fullfile('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/MNI PD25 histology (Xiao 2017) - imported from PD25_MNI',ROI_name); 
    
   %Neudorfer Atlas    
   elseif isequal(Atlas_name, 'Neudorfer')
    
        %Atlas directory
        ROI_dir = fullfile('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/Atlas of the Human Hypothalamus (Neudorfer & Germann 2020)',ROI_name); 
   
   %Custom 
    elseif isequal(Atlas_name, 'Custom')
        
        %Atlas directory 
        ROI_dir = fullfile('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/Custom Drawn ROIs',ROI_name);
        
    end    

    %% Cut voxels to specified space as specified by ROI or grid space

    %Define method for establishing voxel values
    % 1 = consider only the single center voxel BOLD signal
    % 2 = average signal among all voxels at center or adjacent
    voxel_method = 1;

    %Display method specified 
    if isequal(voxel_method, 1)

        %Note: this is equivalent to using all voxels in an ROI
        disp('Method: Single voxel signal method')

    elseif isequal(voxel_method, 2)

        disp('Method: Average adjacent voxels signal method')

    end
    
    %Initialize voxel location variable
    voxelInfo = [];
    xyz = [];
        
    %Load voxelInfo for that ROI or masked
    cd(ROI_dir)
    load(fullfile(ROI_dir,'voxelInfo.mat'));

    %Setup voxel info
    if not(isempty(voxelInfo))

        %Define voxel coordinates
        vXYZ = voxelInfo;

    elseif not(isempty(xyz))

        %Define voxel coordinates     
        vXYZ = xyz.vXYZ; 

    end

    %Specify image size (should be 91x109x91);
    sizeMat = [91,109,91]; 

    %Get the list of index of useful voxels sizeMat(2:4), is the range of X, Y, Z indices. sizeMat(1) is time
    ROI_outInd = mat2ind(vXYZ,sizeMat);  %mat2ind creates an 1D index for the given matrix size (3D)

    %% Load Kmeans ROIs

    %Two versions of the kmeans ROIs - sig voxels only for PP/PnP data
    %and all voxels for CP/CnP data

    %ROI cell
    Kclust_ROI_cell = {'kcluster 1','kcluster 2','kcluster 3'};

    %Loop over clusters
    for kclust = 1:length(Kclust_ROI_cell)

        %ROI name
        kclust_ROI_name = Kclust_ROI_cell{kclust};

        %K-means directory - no-report
        ROI_sig_dir = fullfile('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/CP CnP Kmeans 3 Clusters',[kclust_ROI_name,' noreport voxels']);

        %Load voxelInfo for that ROI or masked
        cd(ROI_sig_dir)
        load(fullfile(ROI_sig_dir,'voxelInfo.mat'));

        %Setup voxel info
        if not(isempty(voxelInfo))

            %Define voxel coordinates
            vXYZ = voxelInfo;

        elseif not(isempty(xyz))

            %Define voxel coordinates     
            vXYZ = xyz.vXYZ; 

        end

        %Specify image size (should be 91x109x91);
        sizeMat = [91,109,91]; 

        %Get the list of index of useful voxels sizeMat(2:4), is the range of X, Y, Z indices. sizeMat(1) is time
        kclust_outInd = mat2ind(vXYZ,sizeMat);  %mat2ind creates an 1D index for the given matrix size (3D)       

        %Rename variable - sig voxels used for PP/PnP data
        eval(['kclust_sigvoxels_',num2str(kclust),'_outInd = kclust_outInd;'])

        %K-means directory - report
        ROI_all_dir = fullfile('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/CP CnP Kmeans 3 Clusters',[kclust_ROI_name,' report voxels']);     

        %Load voxelInfo for that ROI or masked
        cd(ROI_all_dir)
        load(fullfile(ROI_all_dir,'voxelInfo.mat'));

        %Setup voxel info
        if not(isempty(voxelInfo))

            %Define voxel coordinates
            vXYZ = voxelInfo;

        elseif not(isempty(xyz))

            %Define voxel coordinates     
            vXYZ = xyz.vXYZ; 

        end

        %Specify image size (should be 91x109x91);
        sizeMat = [91,109,91]; 

        %Get the list of index of useful voxels sizeMat(2:4), is the range of X, Y, Z indices. sizeMat(1) is time
        kclust_outInd = mat2ind(vXYZ,sizeMat);  %mat2ind creates an 1D index for the given matrix size (3D)       

        %Rename variable - sig voxels used for PP/PnP data
        eval(['kclust_allvoxels_',num2str(kclust),'_outInd = kclust_outInd;'])

    end
        
    %% Load Group BOLD PP PnP data
    
    disp('Loading PP and PnP data')
    
    cd(irrel_group_data_dir)
    load(['Group_PC_BOLD_irrelevant_PP_PnP_',location_data_name,'_threshold_score_thres_',num2str(confidence_score),'_data.mat'])
    
    %Array of data types
    data_type_array = {'group_PP_BOLD_PC_data','group_PnP_BOLD_PC_data'};
    
    %Loop over data types (e.g., main data and subtraction data)
    for type = 1:size(data_type_array, 2)

        %Define current data type
        current_data_type = data_type_array{type};

         %Loop over k-means clusters
        for kclust = 1:length(Kclust_ROI_cell)
            
            %Define k cluster
            kclust_ROI_outInd = eval(['kclust_sigvoxels_',num2str(kclust),'_outInd']);
            
            %Initialize NaN matrix to store BOLD voxel values
            %group_pc_voxel_subset_data = nan(size(ROI_outInd,1),size(eval(current_data_type),1),size(eval(current_data_type),3));
            group_pc_voxel_subset_data = nan(sum(ismember(kclust_ROI_outInd,ROI_outInd)),size(eval(current_data_type),1),size(eval(current_data_type),3));

            %This for loop cycles through all of the participants
            for sess = 1:size(eval(current_data_type),3)

                %This for loop cycles through all time points per participant
                for timepoint = 1:size(eval(current_data_type),1) 

                    %Reshape session voxel data to 91, 109, 91 for each time point and subject
                    img = reshape(eval([current_data_type,'(timepoint,:,sess)']),[91, 109, 91]);

                    %METHOD 1: use only single center voxel value
                    if isequal(voxel_method, 1)

                        %Select subset of voxels [voxels x time x subjects] - Grab the data from all the voxels specified in outInd form the 3D image
                        group_pc_voxel_subset_data(:,timepoint,sess) = img(kclust_ROI_outInd(ismember(kclust_ROI_outInd,ROI_outInd)));%img(ROI_outInd);

                    %METHOD 2: average signal from adjacent voxels
                    elseif isequal(voxel_method, 2)

                        %Loop over voxels
                        for voxel = 1:size(vXYZ,2)

                            %Find the mean BOLD signal for all voxels within spatial
                            %cluster, includes center and adjacent voxels - nanmean([adjacent voxels + center voxel])
                            group_pc_voxel_subset_data(voxel,timepoint,sess) = nanmean([img(all_adjacent_voxels_cell{voxel}); img(ROI_outInd(voxel))]);

                        end

                    end

                end

            end

        %Rename variable to original name
        eval([current_data_type,'_kclust',num2str(kclust),'_ROI = group_pc_voxel_subset_data;']) % this is the ROI percent timecourse averaged across sessions

        end
        
    end
    
    %Clear all voxel variables
    clearvars group_PP_BOLD_PC_data group_PnP_BOLD_PC_data

    %% Load group BOLD CP and CnP data
    
    disp('Loading CP and CnP data')

	cd(rel_group_data_dir)
    load('Group_PC_BOLD_CP_CnP_threshold_data.mat')
    
    %Array of data types
    data_type_array = {'group_CP_BOLD_PC_data','group_CnP_BOLD_PC_data'};
    
    %Loop over data types (e.g., main data and subtraction data)
    for type = 1:size(data_type_array, 2)

        %Define current data type
        current_data_type = data_type_array{type};
        
        %Loop over k-means clusters
        for kclust = 1:length(Kclust_ROI_cell)
            
            %Define k cluster
            kclust_ROI_outInd = eval(['kclust_allvoxels_',num2str(kclust),'_outInd']);
            
            %Initialize NaN matrix to store BOLD voxel values
            %group_pc_voxel_subset_data = nan(size(ROI_outInd,1),size(eval(current_data_type),1),size(eval(current_data_type),3));
            group_pc_voxel_subset_data = nan(sum(ismember(kclust_ROI_outInd,ROI_outInd)),size(eval(current_data_type),1),size(eval(current_data_type),3));

            %This for loop cycles through all of the participants
            for sess = 1:size(eval(current_data_type),3)

                %This for loop cycles through all time points per participant
                for timepoint = 1:size(eval(current_data_type),1) 

                    %Reshape session voxel data to 91, 109, 91 for each time point and subject
                    img = reshape(eval([current_data_type,'(timepoint,:,sess)']),[91, 109, 91]);

                    %METHOD 1: use only single center voxel value
                    if isequal(voxel_method, 1)

                        %Select subset of voxels [voxels x time x subjects] - Grab the data from all the voxels specified in outInd form the 3D image
                        group_pc_voxel_subset_data(:,timepoint,sess) = img(kclust_ROI_outInd(ismember(kclust_ROI_outInd,ROI_outInd)));%img(ROI_outInd);

                    %METHOD 2: average signal from adjacent voxels
                    elseif isequal(voxel_method, 2)

                        %Loop over voxels
                        for voxel = 1:size(vXYZ,2)

                            %Find the mean BOLD signal for all voxels within spatial
                            %cluster, includes center and adjacent voxels - nanmean([adjacent voxels + center voxel])
                            group_pc_voxel_subset_data(voxel,timepoint,sess) = nanmean([img(all_adjacent_voxels_cell{voxel}); img(ROI_outInd(voxel))]);

                        end

                    end

                end

            end

            %Rename variable to original name
            eval([current_data_type,'_kclust',num2str(kclust),'_ROI = group_pc_voxel_subset_data;']) % this is the ROI percent timecourse averaged across sessions

        end
        
    end
    
    %Clear all voxel variables
    clearvars group_CP_BOLD_PC_data group_CnP_BOLD_PC_data
    
    %% CP/CnP PP/PnP and average over voxels and Plot and Run Permutation test
    
    %Number of subjects (PP and PnP Subject number are equal for confidence
    %score threshold of 0.75)
    PP_num_subjects = size(PP_epochs_subjects_list,1);
    CP_num_subjects = size(all_included_subjects_list,1);

    %Time definition 
    timevector = [-20:20];

    %Setup plot
    figure
    hold on

    %Figure parameters
    title(['BOLD ROI Timecourses - ',Atlas_name,' ',ROI_name],'Interpreter','none')
    xlabel('Time (s)')
    ylabel('Percent Change BOLD Signal')
    ylim([-0.15 0.2])
    xlim([-20 20])

    %Plot reference lines
    plot(timevector, zeros(length(timevector)), 'k');
    plot([0 0],[-0.2 0.3], 'k');
    plot([6 6],[-0.2 0.3], '--k'); %Prestim question
    plot([10 10],[-0.2 0.3], '--k'); %Prestim question
    plot([-6 -6],[-0.2 0.3], '--k'); %Poststim question
    plot([-10 -10],[-0.2 0.3], '--k'); %Poststim question
    
    %Specify baseline period (defined between 1 and 41; stimulus onset at 21)
    baseline_window = 1:20;

    %Figure color options
    report_color_cell = {'r','y','b'};
    noreport_color_cell = {'m','g','c'};
    
    %Sig times location on y-axis
    sig_time_cell = [0.19, 0.18, 0.17];
    
    %Subtract main data from subtract data
    
    %Remove subjects without PP and PnP trial instances (subject 603 does
    %not have PP trials)
    if length(PnP_epochs_subjects_list) > length(PP_epochs_subjects_list)
        
        %Cut subject 603 from PnP
        group_PnP_BOLD_PC_data_kclust1_ROI(:,:,not(ismember(PnP_epochs_subjects_list,PP_epochs_subjects_list))) = [];
        group_PnP_BOLD_PC_data_kclust2_ROI(:,:,not(ismember(PnP_epochs_subjects_list,PP_epochs_subjects_list))) = [];
        group_PnP_BOLD_PC_data_kclust3_ROI(:,:,not(ismember(PnP_epochs_subjects_list,PP_epochs_subjects_list))) = [];
       
        %Save list of the subjects tested
        tested_subjects_list = PnP_epochs_subjects_list(ismember(PnP_epochs_subjects_list,PP_epochs_subjects_list));
    
    end
    
    %Loop over clusters
    for kclust = 1:length(Kclust_ROI_cell)
        
        group_PP_PnP_pc_sub_data = eval(['group_PP_BOLD_PC_data_kclust',num2str(kclust),'_ROI']) - eval(['group_PnP_BOLD_PC_data_kclust',num2str(kclust),'_ROI']);
        group_CP_CnP_pc_sub_data = eval(['group_CP_BOLD_PC_data_kclust',num2str(kclust),'_ROI']) - eval(['group_CnP_BOLD_PC_data_kclust',num2str(kclust),'_ROI']);

        %Current matrix (voxels x time x subjects) to (time x subject) by averaging
        %over voxel dimension; Note: that these dimensions are updated from the
        %original data that is time x voxels x subjects
        voxel_avg_group_pc_PP_data = squeeze(nanmean(eval(['group_PP_BOLD_PC_data_kclust',num2str(kclust),'_ROI']),1));
        voxel_avg_group_pc_PnP_data = squeeze(nanmean(eval(['group_PnP_BOLD_PC_data_kclust',num2str(kclust),'_ROI']),1));
        voxel_avg_group_pc_PP_PnP_sub_data = squeeze(nanmean(group_PP_PnP_pc_sub_data,1));

        voxel_avg_group_pc_CP_data = squeeze(nanmean(eval(['group_CP_BOLD_PC_data_kclust',num2str(kclust),'_ROI']),1));
        voxel_avg_group_pc_CnP_data = squeeze(nanmean(eval(['group_CnP_BOLD_PC_data_kclust',num2str(kclust),'_ROI']),1));
        voxel_avg_group_pc_CP_CnP_sub_data = squeeze(nanmean(group_CP_CnP_pc_sub_data,1));

        %Average baseline period from CP - CnP/PP - PnP dataset for each subject
        voxel_avg_sub_PP_PnP_baseline_data = repmat(nanmean(voxel_avg_group_pc_PP_PnP_sub_data(baseline_window,:),1), size(voxel_avg_group_pc_PP_PnP_sub_data,1),1);
        voxel_avg_sub_CP_CnP_baseline_data = repmat(nanmean(voxel_avg_group_pc_CP_CnP_sub_data(baseline_window,:),1), size(voxel_avg_group_pc_CP_CnP_sub_data,1),1);
        
        %% CP vs CnP Testing

        disp('Running Permutation CP/CnP vs PP/PnP Test')

        [clusters, pval, t_sums, permutation_distribution] = permutest_TimeCourses(voxel_avg_group_pc_CP_CnP_sub_data(:,ismember(all_included_subjects_list,PP_epochs_subjects_list)), voxel_avg_group_pc_PP_PnP_sub_data, dependent_samples, ...
            p_threshold, num_permutations, two_sided);

        %Find significant clusters pvalue < 0.05
        sig_clust = find(pval < 0.05);

        %Find the significant time points
        sig_CP_CnP_vs_PP_PnP_time_pts = sort([clusters{sig_clust}]);

        %Save output
        cd(fig_folder)
        save(['pc_timecourse_cluster_CP_CnP_vs_PP_PnP_kclust',num2str(kclust),'.mat'],'clusters','pval','t_sums','permutation_distribution');

        %% CP vs PP Testing

        disp('Running Permutation CP vs PP Test')

        [clusters, pval, t_sums, permutation_distribution] = permutest_TimeCourses(voxel_avg_group_pc_CP_data(:,ismember(all_included_subjects_list,PP_epochs_subjects_list)), voxel_avg_group_pc_PP_data, dependent_samples, ...
            p_threshold, num_permutations, two_sided);

        %Find significant clusters pvalue < 0.05
        sig_clust = find(pval < 0.05);

        %Find the significant time points
        sig_CP_vs_PP_time_pts = sort([clusters{sig_clust}]);

        %Save output
        cd(fig_folder)
        save(['pc_timecourse_cluster_CP_vs_PP_kclust',num2str(kclust),'.mat'],'clusters','pval','t_sums','permutation_distribution');

        %% CnP vs PnP Testing

        disp('Running Permutation CnP vs PnP Test')

        [clusters, pval, t_sums, permutation_distribution] = permutest_TimeCourses(voxel_avg_group_pc_CnP_data(:,ismember(all_included_subjects_list,PP_epochs_subjects_list)), voxel_avg_group_pc_PnP_data, dependent_samples, ...
            p_threshold, num_permutations, two_sided);

        %Find significant clusters pvalue < 0.05
        sig_clust = find(pval < 0.05);

        %Find the significant time points
        sig_CnP_vs_PnP_time_pts = sort([clusters{sig_clust}]);

        %Save output
        cd(fig_folder)
        save(['pc_timecourse_cluster_CnP_vs_PnP_kclust',num2str(kclust),'.mat'],'clusters','pval','t_sums','permutation_distribution');
    
        %% Plot timecourses
        
        %Define current color
        report_current_color = report_color_cell{kclust};
        noreport_current_color = noreport_color_cell{kclust};

        %Plot subject average CP minus CnP (NOTE: using the PP subject number for calculating the SEM)
        CP_minus_CnP = plot(timevector, nanmean(voxel_avg_group_pc_CP_CnP_sub_data,2), report_current_color,'DisplayName',['CP_minus_CnP_kclust_',num2str(kclust)]);
        plot(timevector,nanmean(voxel_avg_group_pc_CP_CnP_sub_data,2)+ (nanstd(voxel_avg_group_pc_CP_CnP_sub_data,0,2)/sqrt(CP_num_subjects)),['--',report_current_color]); 
        plot(timevector,nanmean(voxel_avg_group_pc_CP_CnP_sub_data,2)- (nanstd(voxel_avg_group_pc_CP_CnP_sub_data,0,2)/sqrt(CP_num_subjects)),['--',report_current_color]);

        %Plot subject average PP minus PnP
        PP_minus_PnP = plot(timevector, nanmean(voxel_avg_group_pc_PP_PnP_sub_data,2), noreport_current_color);
        plot(timevector,nanmean(voxel_avg_group_pc_PP_PnP_sub_data,2)+ (nanstd(voxel_avg_group_pc_PP_PnP_sub_data,0,2)/sqrt(PP_num_subjects)),['--',noreport_current_color]); 
        plot(timevector,nanmean(voxel_avg_group_pc_PP_PnP_sub_data,2)- (nanstd(voxel_avg_group_pc_PP_PnP_sub_data,0,2)/sqrt(PP_num_subjects)),['--',noreport_current_color]);

%         %Plot subject average CP
%         CP = plot(timevector, nanmean(voxel_avg_group_pc_CP_data,2), 'b');
%         plot(timevector,nanmean(voxel_avg_group_pc_CP_data,2)+ (nanstd(voxel_avg_group_pc_CP_data,0,2)/sqrt(CP_num_subjects)),'--b'); 
%         plot(timevector,nanmean(voxel_avg_group_pc_CP_data,2)- (nanstd(voxel_avg_group_pc_CP_data,0,2)/sqrt(CP_num_subjects)),'--b');
% 
%         %Plot subject average PP
%         PP = plot(timevector, nanmean(voxel_avg_group_pc_PP_data,2), 'c');
%         plot(timevector,nanmean(voxel_avg_group_pc_PP_data,2)+ (nanstd(voxel_avg_group_pc_PP_data,0,2)/sqrt(PP_num_subjects)),'--c'); 
%         plot(timevector,nanmean(voxel_avg_group_pc_PP_data,2)- (nanstd(voxel_avg_group_pc_PP_data,0,2)/sqrt(PP_num_subjects)),'--c');
% 
%         %Plot subject average CnP
%         CnP = plot(timevector, nanmean(voxel_avg_group_pc_CnP_data,2), 'r');
%         plot(timevector,nanmean(voxel_avg_group_pc_CnP_data,2)+ (nanstd(voxel_avg_group_pc_CnP_data,0,2)/sqrt(CP_num_subjects)),'--r'); 
%         plot(timevector,nanmean(voxel_avg_group_pc_CnP_data,2)- (nanstd(voxel_avg_group_pc_CnP_data,0,2)/sqrt(CP_num_subjects)),'--r');
% 
%         %Plot subject average PnP
%         PnP = plot(timevector, nanmean(voxel_avg_group_pc_PnP_data,2), 'm');
%         plot(timevector,nanmean(voxel_avg_group_pc_PnP_data,2)+ (nanstd(voxel_avg_group_pc_PnP_data,0,2)/sqrt(PP_num_subjects)),'--m'); 
%         plot(timevector,nanmean(voxel_avg_group_pc_PnP_data,2)- (nanstd(voxel_avg_group_pc_PnP_data,0,2)/sqrt(PP_num_subjects)),'--m');

        %Plot significant points
        significant_CP_CnP_vs_PP_PnP = scatter(timevector(sig_CP_CnP_vs_PP_PnP_time_pts), ones(length(sig_CP_CnP_vs_PP_PnP_time_pts),1)*sig_time_cell(kclust), 20, 'filled', report_current_color);
%         significant_CP_vs_PP = scatter(timevector(sig_CP_vs_PP_time_pts), ones(length(sig_CP_vs_PP_time_pts),1)*0.18, 20, 'filled', 'b');
%         significant_CnP_vs_PnP = scatter(timevector(sig_CnP_vs_PnP_time_pts), ones(length(sig_CnP_vs_PnP_time_pts),1)*0.17, 20, 'filled', 'r');

    %legend([CP_minus_CnP, PP_minus_PnP, significant_CP_CnP_vs_PP_PnP], 'CP minus CnP','PP minus PnP','CP-CnP vs PP-PnP Sig', 'Location', 'Southwest')

    %Rename plot headers
    eval(['CP_minus_CnP_kclust_',num2str(kclust),' = CP_minus_CnP;'])
    eval(['PP_minus_PnP_kclust_',num2str(kclust),' = PP_minus_PnP;'])
    eval(['significant_CP_CnP_vs_PP_PnP_kclust_',num2str(kclust),' = significant_CP_CnP_vs_PP_PnP;'])

    end
    
    %Legend
    legend([CP_minus_CnP_kclust_1, PP_minus_PnP_kclust_1, significant_CP_CnP_vs_PP_PnP_kclust_1,CP_minus_CnP_kclust_2, PP_minus_PnP_kclust_2, significant_CP_CnP_vs_PP_PnP_kclust_2,CP_minus_CnP_kclust_3, PP_minus_PnP_kclust_3, significant_CP_CnP_vs_PP_PnP_kclust_3],...
        'CP minus CnP DAS','PP minus PnP DAS','CP-CnP vs PP-PnP Sig DAS','CP minus CnP TPN','PP minus PnP TPN','CP-CnP vs PP-PnP Sig TPN','CP minus CnP DMN','PP minus PnP DMN','CP-CnP vs PP-PnP Sig DMN', 'Location', 'Southwest')
    
    %Save figure
    cd(fig_folder)
    savefig([ROI_name,'_PC_BOLD_timecourses_stats.fig'])    
    
    close

    %Clear unnecessary variables
    clearvars group*
    
end
    