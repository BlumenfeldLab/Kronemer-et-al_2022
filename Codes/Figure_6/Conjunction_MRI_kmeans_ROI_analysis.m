%% Find spatial conjunction between anatomical/functional ROI and kmeans cluster

%This code will load the voxel info for the the ROI and kmeans cluster, and
%compare the voxelInfo information and select the conjunction voxels.

%Written by: Sharif Kronemer
%Date: 7/2/2021
%Modified: 7/20/2021

clear

addpath('/mnt/Data5/Sharif Data F/HNCT fMRI Study/spm12_radiological')
addpath(genpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Extract voxel data/Supplementary functions'));
addpath('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/marsbar-0.44');

%% Prompts

%Select atlas
prompt2 = 'Atlas  [MarsBaR, Stanford, Xiao, Neudorfer, Morel, Custom]: ';
Atlas_name = input(prompt2, 's');

%Select ROIs
if isequal(Atlas_name, 'Stanford')
    
    %ROI cell
    ROI_cell = {'motor_cortex_left_Sensorimotor_01'};%'frontal_gyrus_inferior_parietal_LECN_RECN','middle_frontal_frontal_gyrus_inferior_parietal_LECN_RECN','precuneus_medial_PFC_dorsal_DMN_01_04',...
        %'frontal_cortex_operculum_LECN_01_RECN_01_Visuospatial_03_07','inferior_parietal_gyrus_Visuospatial_02_06', 'inferior_parietal_gyrus_ECN_Visuospatial',...
        %'precuneus_dorsal_DMN_04','thalamus_caudate_putamen_Basal_Ganglia_01_02','cerebellum_LECN_05_RECN_05','frontal_gyrus_LECN_01_RECN_01',...
        %'inferior_parietal_gyrus_LECN_03_RECN_03','middle_frontal_gyrus_LECN_02_RECN_02','SMA_AC_anterior_Salience_03','V1_prim_Visual_01',...
        %'operculum_FEF_Visuospatial_03_07','insula_anterior_Salience_02_05','medial_PFC_dorsal_DMN_01','anterior_Salience','Basal_Ganglia',...
        %'dorsal_DMN','ventral_DMN','Sensorimotor','high_Visual','post_Salience','prim_Visual','Visuospatial','LECN','RECN','right_left_ECN','prim_high_Visual','ant_post_Salience'};

elseif isequal(Atlas_name,'MarsBaR')
    
    %ROI cell
    ROI_cell = {'MNI_Parietal_Sup_L+R_roi'};%{'MNI_Temporal_Mid_R+L_roi','MNI_Hippocampus_ParaHippo_R+L_roi','MNI_Lateral_Occipital_L+R_roi',...
        %'MNI_Temporal_Pole_Mid_Med_L+R_roi','MNI_Temporal_L+R_roi_roi','MNI_Frontal_Sup_Orb_R+L_roi',...
        %'MNI_Cerebelum_Crus1_Crus2_L+R_roi','MNI_Cingulate_L+R','Midbrain_Reticular_Formation_roi',...
        %'Intralaminar_Thalamic_Nuclei_v5_symetric_roi','MNI_Insula_L+R_roi','Thalamus L+R',...
        %'SMA L+R','DMNv2 Region 2','MNI_Cerebelum_Crus1_L+R_roi','MNI_Cerebelum_Crus2_L+R_roi','MNI_Frontal_Mid_L+R_roi',...
        %'MNI_Parietal_Inf_Sup_L+R','MNI_Striatum_L+R_roi','MNI_Vermis_Whole_roi','MNI_Amygdala_L+R_roi','MNI_Frontal_Cortex_L+R',...
        %'MNI_Calcarine_L+R_roi','MNI_Fusiform_L+R_roi'};
    
elseif isequal(Atlas_name,'Xiao')
    
    %ROI cell
    ROI_cell = {'Claustrum','Subthalamus'};
    
elseif isequal(Atlas_name,'Neudorfer')
    
    %ROI cell
    ROI_cell = {'Nucleus Basalis'};
    
elseif isequal(Atlas_name,'Morel')
    
    %ROI cell
    ROI_cell = {'Thalamus body'};
        
elseif isequal(Atlas_name,'Custom')
    
    %ROI cell
    ROI_cell = {'NB_custom_ROI_L+R_roi'};
    
end

%% Directories

%K-means directory
kclust1_report_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/CP CnP Kmeans 3 Clusters/kcluster 1 report voxels';%fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/MRI Analysis/K means Clusters/PP_PnP_trials',...
    %'Whole Brain/Combined Task CP CnP Kclust/PP_PnP_cent_quad_score_0.75_5000/Pos Neg Sig Voxels/Cluster num 3/kcluster 1');
kclust1_noreport_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/CP CnP Kmeans 3 Clusters/kcluster 1 noreport voxels';

kclust2_report_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/CP CnP Kmeans 3 Clusters/kcluster 2 report voxels';%fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/MRI Analysis/K means Clusters/PP_PnP_trials',...
    %'Whole Brain/Combined Task CP CnP Kclust/PP_PnP_cent_quad_score_0.75_5000/Pos Neg Sig Voxels/Cluster num 3/kcluster 2');
kclust2_noreport_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/CP CnP Kmeans 3 Clusters/kcluster 2 noreport voxels';

kclust3_report_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/CP CnP Kmeans 3 Clusters/kcluster 3 report voxels';%fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/MRI Analysis/K means Clusters/PP_PnP_trials',...
    %'Whole Brain/Combined Task CP CnP Kclust/PP_PnP_cent_quad_score_0.75_5000/Pos Neg Sig Voxels/Cluster num 3/kcluster 3');
kclust3_noreport_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/CP CnP Kmeans 3 Clusters/kcluster 3 noreport voxels';

%Save directory
root_save_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/MRI Analysis/Conjunction Analysis/ROI Kmeans Conjunction';
mkdir(root_save_dir)

%% Conjunction Analysis

%REPORT Clusters

%Load kmeans voxel info
load(fullfile(kclust1_report_dir,'voxelInfo.mat'));
kclust_1_report_voxelInfo = voxelInfo;

clearvars voxelInfo

load(fullfile(kclust2_report_dir,'voxelInfo.mat'));
kclust_2_report_voxelInfo = voxelInfo;

clearvars voxelInfo

load(fullfile(kclust3_report_dir,'voxelInfo.mat'));
kclust_3_report_voxelInfo = voxelInfo;

clearvars voxelInfo

%NOREPORT Clusters

%Load kmeans voxel info
load(fullfile(kclust1_noreport_dir,'voxelInfo.mat'));
kclust_1_noreport_voxelInfo = voxelInfo;

clearvars voxelInfo

load(fullfile(kclust2_noreport_dir,'voxelInfo.mat'));
kclust_2_noreport_voxelInfo = voxelInfo;

clearvars voxelInfo

load(fullfile(kclust3_noreport_dir,'voxelInfo.mat'));
kclust_3_noreport_voxelInfo = voxelInfo;

clearvars voxelInfo

% Loop over ROI cell
for roi = 1:length(ROI_cell)
    
    %Select current ROI
    ROI_name = ROI_cell{roi};
    
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
     
   %Morel Atlas
   elseif isequal(Atlas_name, 'Morel')
        
        ROI_dir = fullfile('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/Morel Atlas',ROI_name);
                     
    %Custom 
    elseif isequal(Atlas_name, 'Custom')
        
        %Atlas directory 
        ROI_dir = fullfile('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs',ROI_name);
        
    end  
   
    %Load ROI voxelInfo
    cd(ROI_dir)
    load(fullfile(ROI_dir,'voxelInfo.mat'))
    
    try
        
        ROI_voxelInfo = voxelInfo;
    
    catch
        
        ROI_voxelInfo = xyz.vXYZ;
        
    end
    
    clearvars voxelInfo xyz
    
    %Specify image size (should be 91x109x91);
    sizeMat = [91,109,91]; 

    %Get the list of index of useful voxels sizeMat(2:4), is the range of X, Y, Z indices. sizeMat(1) is time
    ROI_outInd = mat2ind(ROI_voxelInfo,sizeMat);  %mat2ind creates an 1D index for the given matrix size (3D)
    k1_report_outInd = mat2ind(kclust_1_report_voxelInfo,sizeMat);
    k2_report_outInd = mat2ind(kclust_2_report_voxelInfo,sizeMat);
    k3_report_outInd = mat2ind(kclust_3_report_voxelInfo,sizeMat);

    k1_noreport_outInd = mat2ind(kclust_1_noreport_voxelInfo,sizeMat);
    k2_noreport_outInd = mat2ind(kclust_2_noreport_voxelInfo,sizeMat);
    k3_noreport_outInd = mat2ind(kclust_3_noreport_voxelInfo,sizeMat);
    
    %% Conjunction Analysis - Shared voxels between k-means cluster and ROI - REPORT
    
    %Conjunction Voxels - These are the voxels shared between the ROI
    %voxels and those voxels from the CP/CnP k=3 clusters
    ROI_k1_report_conjunction = ismember(ROI_outInd,k1_report_outInd);
    ROI_k2_report_conjunction = ismember(ROI_outInd,k2_report_outInd);
    ROI_k3_report_conjunction = ismember(ROI_outInd,k3_report_outInd);
    ROI_unassigned_conjunction = ~ismember(ROI_outInd,k1_report_outInd)& ~ismember(ROI_outInd,k2_report_outInd)& ~ismember(ROI_outInd,k3_report_outInd);
    
    %Conjunction Voxel locations
    ROI_k1_report_voxelinfo = ROI_voxelInfo(:,ROI_k1_report_conjunction);
    ROI_k2_report_voxelinfo = ROI_voxelInfo(:,ROI_k2_report_conjunction);
    ROI_k3_report_voxelinfo = ROI_voxelInfo(:,ROI_k3_report_conjunction);
    ROI_unassigned_voxelinfo = ROI_voxelInfo(:,ROI_unassigned_conjunction);

    %Number of voxel/percent of voxel owned by kmeans clusters
    ROI_k1_report_voxel_num = sum(ismember(ROI_outInd,k1_report_outInd));
    ROI_k1_report_voxel_percent_of_allROIvoxels = (ROI_k1_report_voxel_num/length(ROI_outInd))*100;
    
    ROI_k2_report_voxel_num = sum(ismember(ROI_outInd,k2_report_outInd));
    ROI_k2_report_voxel_percent_of_allROIvoxels = (ROI_k2_report_voxel_num/length(ROI_outInd))*100;
    
    ROI_k3_report_voxel_num = sum(ismember(ROI_outInd,k3_report_outInd));
    ROI_k3_report_voxel_percent_of_allROIvoxels = (ROI_k3_report_voxel_num/length(ROI_outInd))*100;
    
    ROI_unassigned_voxel_num = sum(~ismember(ROI_outInd,k1_report_outInd)& ~ismember(ROI_outInd,k2_report_outInd)& ~ismember(ROI_outInd,k3_report_outInd));
    ROI_unassigned_voxel_percent_of_allROIvoxels = ROI_unassigned_voxel_num/length(ROI_outInd)*100;
    
    %Calculate percent relevant to the number of sig voxels - in other
    %words the total number of voxels in the k means ROIs and not the
    %anatomical or functional ROI which is calculated from the all ROI
    %option above. 
    total_report_kcluster_voxel_num = sum(ROI_k1_report_voxel_num + ROI_k2_report_voxel_num + ROI_k3_report_voxel_num);
    ROI_k1_report_voxel_percent_of_allkclustervoxels = (ROI_k1_report_voxel_num/total_report_kcluster_voxel_num)*100;
    ROI_k2_report_voxel_percent_of_allkclustervoxels = (ROI_k2_report_voxel_num/total_report_kcluster_voxel_num)*100;
    ROI_k3_report_voxel_percent_of_allkclustervoxels = (ROI_k3_report_voxel_num/total_report_kcluster_voxel_num)*100;

   %% Conjunction Analysis - Shared voxels between k-means cluster and ROI - NO-REPORT
    
    %Conjunction Voxels - These are the voxels shared between the ROI
    %voxels and those voxels from the CP/CnP k=3 clusters
    ROI_k1_noreport_conjunction = ismember(ROI_outInd,k1_noreport_outInd);
    ROI_k2_noreport_conjunction = ismember(ROI_outInd,k2_noreport_outInd);
    ROI_k3_noreport_conjunction = ismember(ROI_outInd,k3_noreport_outInd);
    
    %Conjunction Voxel locations
    ROI_k1_noreport_voxelinfo = ROI_voxelInfo(:,ROI_k1_noreport_conjunction);
    ROI_k2_noreport_voxelinfo = ROI_voxelInfo(:,ROI_k2_noreport_conjunction);
    ROI_k3_noreport_voxelinfo = ROI_voxelInfo(:,ROI_k3_noreport_conjunction);

    %Number of voxel/percent of voxel owned by kmeans clusters
    ROI_k1_noreport_voxel_num = sum(ismember(ROI_outInd,k1_noreport_outInd));
    ROI_k1_noreport_voxel_percent_of_allROIvoxels = (ROI_k1_noreport_voxel_num/length(ROI_outInd))*100;
    
    ROI_k2_noreport_voxel_num = sum(ismember(ROI_outInd,k2_noreport_outInd));
    ROI_k2_noreport_voxel_percent_of_allROIvoxels = (ROI_k2_noreport_voxel_num/length(ROI_outInd))*100;
    
    ROI_k3_noreport_voxel_num = sum(ismember(ROI_outInd,k3_noreport_outInd));
    ROI_k3_noreport_voxel_percent_of_allROIvoxels = (ROI_k3_noreport_voxel_num/length(ROI_outInd))*100;
    
    %Calculate percent relevant to the number of sig voxels - in other
    %words the total number of voxels in the k means ROIs and not the
    %anatomical or functional ROI which is calculated from the all ROI
    %option above. 
    total_noreport_kcluster_voxel_num = sum(ROI_k1_noreport_voxel_num + ROI_k2_noreport_voxel_num + ROI_k3_noreport_voxel_num);
    ROI_k1_noreport_voxel_percent_of_allkclustervoxels = (ROI_k1_noreport_voxel_num/total_noreport_kcluster_voxel_num)*100;
    ROI_k2_noreport_voxel_percent_of_allkclustervoxels = (ROI_k2_noreport_voxel_num/total_noreport_kcluster_voxel_num)*100;
    ROI_k3_noreport_voxel_percent_of_allkclustervoxels = (ROI_k3_noreport_voxel_num/total_noreport_kcluster_voxel_num)*100;  
    
    
    %% Create nii files of conjunction
    
    %Kmeans_ROI cells
    kmeans_ROIs = {'ROI_k1_report','ROI_k2_report','ROI_k3_report',...
        'ROI_k1_noreport','ROI_k2_noreport','ROI_k3_noreport','ROI_unassigned'};
    
    %Define ROI transform matrix
    V = transform_matrix;
    
    %Loop over kmeans_ROIs  
    for n = 1:length(kmeans_ROIs)
        
        %Check if voxel info is empty
        if isempty(eval([kmeans_ROIs{n},'_voxelinfo']))
           
            disp(['No voxels for ',kmeans_ROIs{n},' - skipping'])
            continue
            
        end    

        %maroi_pointlist is the Marsbar function to create an ROI
        act_roi = maroi_pointlist(struct('XYZ', eval([kmeans_ROIs{n},'_voxelinfo']), 'mat', V), 'vox'); 

        %Make save directory
        save_dir = fullfile(root_save_dir,Atlas_name,ROI_name);
        mkdir(save_dir)
        
        %Save ROI to MarsBaR ROI file
        cd(save_dir)
        saveroi(act_roi, [kmeans_ROIs{n},'_conjunction.mat']);
        save_as_image(act_roi, [kmeans_ROIs{n},'_conjunction.nii']);
    
    end
    
    %Save voxel stat info
    cd(save_dir)
    save voxel_num_percent_values.mat ROI_k1* ROI_k2* ROI_k3* ROI_unassigned*
    
end
