%% MRI Permutation Conjunction Analysis - CP/CnP vs PP/PnP - Predicted perceived vs true class

%Conjunction test: CP_minus_CnP vs PP_minus_PnP

%Written by: Sharif Kronemer
%Date: 5/27/2021
%Modified: 8/26/2021

clear

%% Prompts

%Run location
run_location = 'l';

%Confidence score
confidence_score = 0.75;

%Only include voxel that are not statistically different between CP-CnP vs
%PP-PnP (1 = yes - non-sig voxels; 0 = no - all voxels)
non_sig_voxels_only = 1;

%Report vs No-Report (1) or No-Report vs Report (2)
comparison_direction = 1; 

%% Directories

addpath('Z:\Sharif Data F\HNCT fMRI Study\spm12_radiological')
addpath(genpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\Extract voxel data\Supplementary functions'));

%CP minus CnP Stats
CP_minus_CnP_dir = 'S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Group Analysis\MRI Analysis\Percent Change Permutation Stats\Whole Brain_CP_minus_CnP_relevant_cent_quad_threshold_perm_5000';

%PP minus PnP Stats
PP_minus_PnP_dir = ['S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Group Analysis\MRI Analysis\Percent Change Permutation Stats\Whole Brain_PP_minus_PnP_irrelevant_cent_quad_threshold_score_thres_',num2str(confidence_score),'_perm_5000'];

if isequal(comparison_direction, 1)

    %CP minus CnP minus PP minus PnP
    CP_CnP_minus_PP_PnP = ['S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Group Analysis\MRI Analysis\Percent Change Permutation Stats\Whole Brain CP_minus_CnP minus PP_minus_PnP perm 5000 con score ',num2str(confidence_score)];

elseif isequal(comparison_direction, 2)
    
    %PP minus PnP minus CP minus CnP
    CP_CnP_minus_PP_PnP = ['S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Group Analysis\MRI Analysis\Percent Change Permutation Stats\Whole Brain PP_minus_PnP minus CP_minus_CnP perm 5000 con score ',num2str(confidence_score)];

end

%Save directory
save_dir = 'S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Group Analysis\MRI Analysis\Conjunction Analysis';

%Grid directory
grid_dir = 'S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\Permutation Statistical Analysis\Brain Voxel Subset\Whole Brain'; 

%MNI template
template_image = 'S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\spm12_radiological\canonical\single_subj_T1.nii';
%template_image = 'R:\RNS Study\RNS Locations\BigBrain Atlas\PD25-MPRAGET1w-template-500um.nii';

%ROI name
if isequal(comparison_direction, 1)
    
    if isequal(non_sig_voxels_only, 1)

        ROI_name = ['Whole Brain CP-CnP vs PP-PnP Score ', num2str(confidence_score), ' CP_CnP_minus_PP_PnP_not_sig_diff_voxels'];

    else

        ROI_name = ['Whole Brain CP-CnP vs PP-PnP Score ', num2str(confidence_score), ' all_voxels'];

    end
    
elseif isequal(comparison_direction, 2)
    
    if isequal(non_sig_voxels_only, 1)

        ROI_name = ['Whole Brain PP-PnP vs CP-CnP Score ', num2str(confidence_score), ' PP_PnP_minus_CP_CnP_not_sig_diff_voxels'];

    else

        ROI_name = ['Whole Brain PP-PnP vs CP-CnP Score ', num2str(confidence_score), ' all_voxels'];

    end
    
end

%Figure folders
fig_img_folder = fullfile(save_dir, ROI_name, ['Cluster_maps']);
hdr_img_folder = fullfile(save_dir, ROI_name, ['Source_imgs']);

%Create figure folders
mkdir(fig_img_folder)
mkdir(hdr_img_folder)

%% Prepare Brain Voxel Information 

%Load voxelInfo for that ROI or masked
cd(grid_dir)
load(fullfile(grid_dir,'voxelInfo.mat'));

%Unique filenames for the constume built ROIs 
try
   
   %Define voxel coordinates
   vXYZ = voxelInfo;
   
 
end

try

    %Define voxel coordinates     
    vXYZ = xyz.vXYZ; 

end

%If use grid voxels
disp('Cut voxels to subset!')

%Load grid_dir
cd(grid_dir)

%Load grid voxels
try

    %Whole brain grid voxel subset
    load('brain_grid_voxel_info.mat')

catch

    %ROI grid voxel subset
    load('ROI_grid_voxel_info.mat')

end

%Cut voxels to grid voxels
vXYZ = vXYZ(:, brain_grid_voxels); 

%Specify image size (should be 91x109x91);
sizeMat = [91,109,91]; 

%Get the list of index of useful voxels sizeMat(2:4), is the range of X, Y, Z indices. sizeMat(1) is time
outInd = mat2ind(vXYZ,sizeMat); 

%% Find Statistically Significant Clusters

%Confirmed Perceived minus Confirmed not Perceived (True class)
cd(CP_minus_CnP_dir)
load('pc_cluster_sumt_rand_005.mat', 'clust_info')

%Find significant clusters pvalue < 0.05 - Negative and positive clusters independently
true_pos_clust = find(clust_info.pos_clust_pval < 0.05);
true_neg_clust = find(clust_info.neg_clust_pval < 0.05);

%Create empty matrix channels by time
true_sig_matrix = zeros(41,902629);

%Loop over postive clusters
for clust = 1:length(true_pos_clust)

    %Loop over time
    for time = 1:size(true_sig_matrix,1)

        %Find postive significant cluster number/idx in voxel x time matrix
        [sig_grid_voxels, col] = find(clust_info.pos_clust_ids(:,time) == true_pos_clust(clust));

        %Include voxel adjacency in significant voxel designation (output = a cell of adjacent voxels for each of the significant voxels)
        adj_voxels = all_adjacent_voxels_cell(sig_grid_voxels);

        %Initialize array variable used to covert adjacent cell to array
        adj_voxels_array = [];

        %Loop over sig_grid_voxels
        for voxel = 1:size(sig_grid_voxels,1)

            %Store all adjacent voxels for each sig_grid_voxels in array
            adj_voxels_array = [adj_voxels_array; cell2mat(adj_voxels(voxel))];

        end

        %Apply constant for plotting.

        %Give post clusters constant value - find the voxel all volume idx
        true_sig_matrix(time,outInd(sig_grid_voxels)) = 1;

        %Give adjacent voxels to post clusters constant value
        true_sig_matrix(time,adj_voxels_array) = 1;      

    end

end

%Loop over negative clusters
for clust = 1:length(true_neg_clust)

    %Loop over time
    for time = 1:size(true_sig_matrix,1)

        %Find postive significant cluster number/idx in voxel x time matrix
        [sig_grid_voxels, col] = find(clust_info.neg_clust_ids(:,time) == true_neg_clust(clust));

        %Include voxel adjacency in significant voxel designation (output = a cell of adjacent voxels for each of the significant voxels)
        adj_voxels = all_adjacent_voxels_cell(sig_grid_voxels);

        %Initialize array variable used to covert adjacent cell to array
        adj_voxels_array = [];

        %Loop over sig_grid_voxels
        for voxel = 1:size(sig_grid_voxels,1)

            %Store all adjacent voxels for each sig_grid_voxels in array
            adj_voxels_array = [adj_voxels_array; cell2mat(adj_voxels(voxel))];

        end

        %Apply constant for plotting.

        %Give neg clusters constant value - find the voxel all volume idx
        true_sig_matrix(time,outInd(sig_grid_voxels)) = -1;

        %Give adjacent voxels to neg clusters constant value
        true_sig_matrix(time,adj_voxels_array) = -1;      

    end

end

%Predicted Perceived minus Predicted not Perceived - Predicted Class
cd(PP_minus_PnP_dir)
load('pc_cluster_sumt_rand_005.mat', 'clust_info')

%Find significant clusters pvalue < 0.05 - Negative and positive clusters independently
predicted_pos_clust = find(clust_info.pos_clust_pval < 0.05);
predicted_neg_clust = find(clust_info.neg_clust_pval < 0.05);

%Create empty matrix channels by time
predicted_sig_matrix = zeros(41,902629);

%Loop over postive clusters
for clust = 1:length(predicted_pos_clust)

    %Loop over time
    for time = 1:size(predicted_sig_matrix,1)

        %Find postive significant cluster number/idx in voxel x time matrix
        [sig_grid_voxels, col] = find(clust_info.pos_clust_ids(:,time) == predicted_pos_clust(clust));

        %Include voxel adjacency in significant voxel designation (output = a cell of adjacent voxels for each of the significant voxels)
        adj_voxels = all_adjacent_voxels_cell(sig_grid_voxels);

        %Initialize array variable used to covert adjacent cell to array
        adj_voxels_array = [];

        %Loop over sig_grid_voxels
        for voxel = 1:size(sig_grid_voxels,1)

            %Store all adjacent voxels for each sig_grid_voxels in array
            adj_voxels_array = [adj_voxels_array; cell2mat(adj_voxels(voxel))];

        end

        %Apply constant for plotting.

        %Give post clusters constant value - find the voxel all volume idx
        predicted_sig_matrix(time,outInd(sig_grid_voxels)) = 1;

        %Give adjacent voxels to post clusters constant value
        predicted_sig_matrix(time,adj_voxels_array) = 1;      

    end

end

%Loop over negative clusters
for clust = 1:length(predicted_neg_clust)

    %Loop over time
    for time = 1:size(predicted_sig_matrix,1)

        %Find postive significant cluster number/idx in voxel x time matrix
        [sig_grid_voxels, col] = find(clust_info.neg_clust_ids(:,time) == predicted_neg_clust(clust));

        %Include voxel adjacency in significant voxel designation (output = a cell of adjacent voxels for each of the significant voxels)
        adj_voxels = all_adjacent_voxels_cell(sig_grid_voxels);

        %Initialize array variable used to covert adjacent cell to array
        adj_voxels_array = [];

        %Loop over sig_grid_voxels
        for voxel = 1:size(sig_grid_voxels,1)

            %Store all adjacent voxels for each sig_grid_voxels in array
            adj_voxels_array = [adj_voxels_array; cell2mat(adj_voxels(voxel))];

        end

        %Apply constant for plotting.

        %Give neg clusters constant value - find the voxel all volume idx
        predicted_sig_matrix(time,outInd(sig_grid_voxels)) = -1;

        %Give adjacent voxels to neg clusters constant value
        predicted_sig_matrix(time,adj_voxels_array) = -1;      

    end

end

%Select only non-sig voxels from CP-CnP vs PP-PnP stats
if isequal(non_sig_voxels_only,1)

    %CP minus CnP vs PP minus PnP
    cd(CP_CnP_minus_PP_PnP)
    load('pc_cluster_sumt_rand_005.mat', 'clust_info')

    %Find significant clusters pvalue < 0.05 - Negative and positive clusters independently
    true_vs_pred_pos_clust = find(clust_info.pos_clust_pval < 0.05);
    true_vs_pred_neg_clust = find(clust_info.neg_clust_pval < 0.05);

    %Create empty matrix channels by time
    true_vs_pred_sig_matrix = zeros(41,902629);

    %Loop over postive clusters
    for clust = 1:length(true_vs_pred_pos_clust)

        %Loop over time
        for time = 1:size(true_vs_pred_sig_matrix,1)

            %Find postive significant cluster number/idx in voxel x time matrix
            [sig_grid_voxels, col] = find(clust_info.pos_clust_ids(:,time) == true_vs_pred_pos_clust(clust));

            %Include voxel adjacency in significant voxel designation (output = a cell of adjacent voxels for each of the significant voxels)
            adj_voxels = all_adjacent_voxels_cell(sig_grid_voxels);

            %Initialize array variable used to covert adjacent cell to array
            adj_voxels_array = [];

            %Loop over sig_grid_voxels
            for voxel = 1:size(sig_grid_voxels,1)

                %Store all adjacent voxels for each sig_grid_voxels in array
                adj_voxels_array = [adj_voxels_array; cell2mat(adj_voxels(voxel))];

            end

            %Apply constant for plotting.

            %Give post clusters constant value - find the voxel all volume idx
            true_vs_pred_sig_matrix(time,outInd(sig_grid_voxels)) = 1;

            %Give adjacent voxels to post clusters constant value
            true_vs_pred_sig_matrix(time,adj_voxels_array) = 1;      

        end

    end

    %Loop over negative clusters
    for clust = 1:length(true_vs_pred_neg_clust)

        %Loop over time
        for time = 1:size(true_vs_pred_sig_matrix,1)

            %Find postive significant cluster number/idx in voxel x time matrix
            [sig_grid_voxels, col] = find(clust_info.neg_clust_ids(:,time) == true_vs_pred_neg_clust(clust));

            %Include voxel adjacency in significant voxel designation (output = a cell of adjacent voxels for each of the significant voxels)
            adj_voxels = all_adjacent_voxels_cell(sig_grid_voxels);

            %Initialize array variable used to covert adjacent cell to array
            adj_voxels_array = [];

            %Loop over sig_grid_voxels
            for voxel = 1:size(sig_grid_voxels,1)

                %Store all adjacent voxels for each sig_grid_voxels in array
                adj_voxels_array = [adj_voxels_array; cell2mat(adj_voxels(voxel))];

            end

            %Apply constant for plotting.

            %Give neg clusters constant value - find the voxel all volume idx
            true_vs_pred_sig_matrix(time,outInd(sig_grid_voxels)) = -1;

            %Give adjacent voxels to neg clusters constant value
            true_vs_pred_sig_matrix(time,adj_voxels_array) = -1;      

        end

    end

end

%% Conjunction Analysis 

%Conjunction matrix
conjunction_matrix = zeros(41,902629);

%Add the true and predicted sig matrixes identifying common sig voxels in
%negative and positive clusters
additive_matrix = true_sig_matrix + predicted_sig_matrix;

%Positive cluster conjunction (value of 2 indicates two postive clusters in
%time and space)

%All voxels
if isequal(non_sig_voxels_only,0)

    conjunction_matrix(additive_matrix == 2) = 1;%3; 

else

    %Only include voxels not sig different between CP-CnP vs PP-PnP
    conjunction_matrix(find(additive_matrix == 2 & true_vs_pred_sig_matrix == 0)) = 1;%3; 

end

%Negative cluster conjunction (value of -2 indicates two negative clusters
%in time and space)

%All voxels
if isequal(non_sig_voxels_only,0)

    conjunction_matrix(additive_matrix == -2) = -1; 

else

    %Only include voxels not sig different between CP-CnP vs PP-PnP
    conjunction_matrix(find(additive_matrix == -2 & true_vs_pred_sig_matrix == 0)) = -1; 

end
    
%% Plot Brain Map of Positive and Negative Conjunction

%Rename variable dimensions - time x voxel
mean_2D_PC_cut = conjunction_matrix;

%Loop over time points 
%data_pos = 1x41 matrix (1xtime) each time point has a cell of the voxels
%in XYZ space
%data_neg = 1x41 matrix (1xtime)  each time point has a cell of the voxels
%in XYZ space
for time = 1:size(mean_2D_PC_cut,1)

    %Reshape data into 91x109x91 matrix from 1x902629
    img_3D_pertime = reshape(mean_2D_PC_cut(time,:),[91, 109, 91]);

    %Initialize variable 91x109x91
    data_pos{time} = zeros(size(img_3D_pertime));
    data_neg{time} = zeros(size(img_3D_pertime));

    %Find positive and negative values
    data_pos{time} = (img_3D_pertime > 0.01).*img_3D_pertime; %Percent change values above 0.01 are consider positive
    data_neg{time} = (img_3D_pertime < -0.01).*img_3D_pertime; %Percent cahnge values below -0.01 are considered negative

end

%Open template voxelInfo for MNI brain
if isequal(run_location, 'l')

    load('Z:\Sharif Data F\HNCT fMRI Study\fMRI Analysis\Subject Analysis\191AM\MRI Analysis\Extracted voxel data\Run 1 Movie\axialSlices\voxelInfo.mat')

elseif isequal(run_location, 's')

    load('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Subject Analysis/191AM/MRI Analysis/Extracted voxel data/Run 1 Movie/axialSlices/voxelInfo.mat')

end

%Define XYZ data for template image
xyz = xyz.vXYZ; %Load axial slices MR data 3x187076
avw = load_nii(template_image); %Canonical MNI space 91x109x91

%Loop over each timepoint (41)
for time = 1:size(data_pos, 2)

    %Convert the current postive data at time point to double
    newimg = double(data_pos{time});

    %Multiple template image data by 0
    avw.img = avw.img*0;

    %Loop over voxels
    for i = 1: size(xyz, 2)

        avw.img(xyz(1,i), xyz(2, i), xyz(3, i)) = newimg(size(newimg,1)-xyz(1,i)+1,xyz(2, i),xyz(3, i)); %flip L-R

    end

    avw.hdr.dime.datatype = 16;
    avw.hdr.dime.bitpix = int16(32);
    avw.hdr.dime.cal_max = max(max(max(avw.img)));
    avw.hdr.dime.roi_scale = 1;

    %Save nii file
    cd(hdr_img_folder)
    save_nii(avw, ['Act_time_' num2str(time)]); 
    
    %Convert the current negative data at time point to double
    newimg = double(data_neg{time});

    %Multiple template image data by 0
    avw.img = avw.img*0;

    %Loop over voxels
    for i = 1: size(xyz, 2)

        avw.img(xyz(1,i), xyz(2, i), xyz(3, i)) = -newimg(size(newimg,3)-xyz(1,i)+1, xyz(2, i), xyz(3, i)); %flip L-R

    end

    avw.hdr.dime.datatype = 16;
    avw.hdr.dime.bitpix = int16(32);
    avw.hdr.dime.cal_max = max(max(max(avw.img)));
    avw.hdr.dime.roi_scale = 1;

    %Save nii file
    cd(hdr_img_folder)
    save_nii(avw, ['Dec_time_' num2str(time)]);

    end

    %Setup SPM Model
    global model;
    model.xacross = 'auto';
    model.itype{1} = 'Structural';
    model.itype{2} = 'Blobs - Positive';
    model.itype{3} = 'Blobs - Negative';
    model.imgns{1} = 'Img 1 (Structural)';
    model.imgns{2} = 'Img 2 (Blobs - Positive)';
    model.imgns{3} = 'Img 3 (Blobs - Negative)';
    model.range(:,1) = [0 1];

    %Define plotting scale 
    model.range(:,2) = [0.15; 3]; 
    model.range(:,3) = [0.15; 3]; 

    model.transform = 'axial'; %'axial','coronal','sagittal'
    model.axialslice = [-48:8:62];%[-56:8:86];%[-26:2:56]; %[-56:6:86] Defines the # of slices and slice thickness plotted; full range:[-72:2:102]
    model.coronalslice = [-92:8:52]; %Defines the # of slices and slice thickness plotted; full range:[-72:2:102]

    %Loop over time points
    for time = 1:size(data_pos, 2)

    %Setup model
    model.imgs{1, 1} = template_image;
    model.imgs{2, 1} = fullfile(hdr_img_folder, ['Act_time_' num2str(time) '.img']);
    model.imgs{3, 1} = fullfile(hdr_img_folder, ['Dec_time_' num2str(time) '.img']);

    %Run image generation
    display_slices_bai_1;

    cd(fig_img_folder);
    set(gcf,'PaperPositionMode','auto');
    print('-dtiff', ['Time_' num2str(time)]);
    cd('..');
    close all;

end
