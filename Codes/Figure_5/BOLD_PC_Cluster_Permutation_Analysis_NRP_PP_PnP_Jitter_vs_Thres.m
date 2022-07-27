%% MR Voxel Based Permutation Analysis - PP and PnP Data - Jitter vs Threshold

%This code will:

%(1) Load group PP and PnP PC BOLD data
%(2) Cut voxels to ROI size (Whole brain voxel subset)
%(3) Subtract PP/PnP data if specified and baseline data to prestimulus
%period (1-20 seconds; stimulus onset time is 21)
%(4) Run cluster based permutation analysis 
%(5) Generate plots of the signifcant voxels on a brain map according to
%the specified visualization method and save image files

%Note: The run location has been hard coded to the server because of
%computional inefficiency 

%Written by: Sharif I. Kronemer
%Date: 6/4/2021
%Last Modified: 6/10/2022

clear all

%% Prompts

% %Run location
% prompt1 = 'Run location (s/l): ';
% run_location = input(prompt1, 's');
run_location = 's'; %Server

%Select data type
prompt3 = 'Threshold stimulus to analyze (PP,PnP,PP_minus_PnP): ';
thres_data_type = input(prompt3, 's');

%Select data type
prompt3 = 'Jitter stimulus to analyze (PP,PnP,PP_minus_PnP): ';
jitter_data_type = input(prompt3, 's');

%Report minus No report or No report minus report 
prompt5 = 'Subtraction direction (No_Report, ISI): ';
sub_direction = input(prompt5, 's');

%% Parameters

%Number of permutations
num_permutation = 5000;

%Confidence Score
confidence_score = 0.75;

%Center and quadrant name variable (e.g., 'cent','quad', and 'cent_quad')
location_name = 'cent_quad';

%Irrelevant/Relevant conditions (relevant/irrelevant/rel_irrel)
task_condition = 'rel_irrel';

%All jitter trials (1) or CnP/matched jitter trials (0)
all_jitter_trials = 0;

%% ROI Directories 

%Name ROI
ROI_string = 'Whole Brain';

%Save directory name
if isequal(all_jitter_trials, 1)
    
    if isequal(sub_direction,'No_Report')
    
        ROI_name = [ROI_string,'_thres_',thres_data_type,'_vs_jitter_',jitter_data_type,'_',task_condition,'_',location_name,'_',num2str(confidence_score),'_perm_',num2str(num_permutation)];

    elseif isequal(sub_direction,'ISI')
        
        ROI_name = [ROI_string,'_jitter_',thres_data_type,'_vs_thres_',jitter_data_type,'_',task_condition,'_',location_name,'_',num2str(confidence_score),'_perm_',num2str(num_permutation)];

    end
    
else
   
    if isequal(sub_direction,'No_Report')
        
        ROI_name = [ROI_string,'_thres_',thres_data_type,'_vs_jitter_matched_',jitter_data_type,'_',task_condition,'_',location_name,'_',num2str(confidence_score),'_perm_',num2str(num_permutation)];

    elseif isequal(sub_direction,'ISI')
        
        ROI_name = [ROI_string,'_jitter_matched_',thres_data_type,'_vs_thres_',jitter_data_type,'_',task_condition,'_',location_name,'_',num2str(confidence_score),'_perm_',num2str(num_permutation)];
        
    end
    
end

%ROI file directory

%Local
if isequal(run_location, 'l')

elseif isequal(run_location, 's')
    
    %ROI_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/Marsbar ROIs/MNI_Parietal_Sup_L_roi'; %Parietal
    %ROI_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/Marsbar ROIs/Intralaminar_Thalamic_Nuclei_v5_symetric_roi'; %Thalamic
    %ROI_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/Marsbar ROIs/Midbrain_Reticular_Formation_roi';
    %ROI_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/FEF ROI Custom/voxelInfo 91 109 91';
    %ROI_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/Marsbar ROIs/MNI_Frontal_Mid_L+R_roi_roi';
    %ROI_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Formatted ROIs/FEF ROI Custom/voxelInfo 91 109 91';
    
    %Whole Brain
    %ROI_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Subject Analysis/191AM/MRI Analysis/Extracted voxel data/Run 1 Movie/axialSlices'; %MNI whole brain
    ROI_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Permutation Statistical Analysis/Brain Voxel Subset', ROI_string); %Any ROI inputted above
       
end

%% Add Functions to Path

%Local
if isequal(run_location, 'l')

%Server
elseif isequal(run_location, 's')
    
    %Add SPM and supplementary functions to path
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Permutation Statistical Analysis')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Permutation Statistical Analysis/Permutation functions')
    addpath('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Time courses/Supplementary functions/fMRI Analysis/Code for Percent Change Maps/subRoutines')
    addpath('/mnt/Data5/Sharif Data F/HNCT fMRI Study/spm12_radiological')
    addpath(genpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Extract voxel data/Supplementary functions'));
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/Behavioral Analysis')
    
    %Brain grid info
    grid_dir = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Permutation Statistical Analysis/Brain Voxel Subset/', ROI_string); %Any ROI inputted above

    %MNI brain
    template_image = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/spm12_radiological/canonical/single_subj_T1.nii';

    %Save diretory
    save_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Analysis/MRI Analysis/Percent Change Permutation Stats';
    
    %Group data directories
    data_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Group Data/MRI/BOLD PC';
    
    %Figure folders
    fig_folder = fullfile(save_dir, ROI_name);

    %Create figure folders
    mkdir(fig_folder)

end

%% Load group PP and PnP data [time x voxel x subjects]

%Load PP and PnP data
cd(data_dir)

%Jitter Data
if isequal(all_jitter_trials, 1)
    
    load(['Group_PC_BOLD_',task_condition,'_PP_PnP_',location_name,'_jitter_score_thres_',num2str(confidence_score),'_data.mat'])

else
    
    load(['Group_PC_BOLD_',task_condition,'_PP_PnP_',location_name,'_jitter_matched_score_thres_',num2str(confidence_score),'_data.mat'])
    
end

%Rename variables
jitter_PP_BOLD_PC_data = group_PP_BOLD_PC_data;
jitter_PnP_BOLD_PC_data = group_PnP_BOLD_PC_data;
jitter_PP_epochs_subjects_list = PP_epochs_subjects_list;
jitter_PnP_epochs_subjects_list = PnP_epochs_subjects_list;

%Clear variables to save memory
clearvars group_PP* group_PnP*

%Threshold Data
load(['Group_PC_BOLD_irrelevant_PP_PnP_',location_name,'_threshold_score_thres_',num2str(confidence_score),'_data.mat'])

%Rename variables
thres_PP_BOLD_PC_data = group_PP_BOLD_PC_data;
thres_PnP_BOLD_PC_data = group_PnP_BOLD_PC_data;
thres_PP_epochs_subjects_list = PP_epochs_subjects_list;
thres_PnP_epochs_subjects_list = PnP_epochs_subjects_list;

%Clear variables to save memory
clearvars group_PP* group_PnP*

%% Cut voxels to specified space as specified by ROI or grid space

%Define method for establishing voxel values
% 1 = consider only the single center voxel BOLD signal
% 2 = average signal among all voxels at center or adjacent
voxel_method = 2;

% %Use grid voxels or all voxels
% voxel_method = 1; %1 = USE grid voxels; 0 = do not use grid, instead use all voxels

%Display method specified 
if isequal(voxel_method, 1) %Use all voxels
    
    disp('Method: Single voxel signal method')
    
    %Voxel proximity threshold distance - used for voxel neighborhood 
    distance_threshold = 1; %Note: When using grid voxels distance is 3, but when considering all voxels distance of 1 is optimal

elseif isequal(voxel_method, 2) %Use grid/subset voxels
    
    disp('Method: Average adjacent voxels signal method')
    
    %Voxel proximity threshold distance - used for voxel neighborhood
    distance_threshold = 3; %Note: When using grid voxels distance is 3, but when considering all voxels distance of 1 is optimal

end

%PP minus PnP compairson
if isequal(jitter_data_type, 'PP_minus_PnP') && isequal(thres_data_type,'PP_minus_PnP')

    %Array of data types
    data_type_array = {'jitter_PP_BOLD_PC_data', 'jitter_PnP_BOLD_PC_data', 'thres_PP_BOLD_PC_data', 'thres_PnP_BOLD_PC_data'};

% PP compairson
elseif isequal(jitter_data_type, 'PP') && isequal(thres_data_type,'PP')

    %Array of data types
    data_type_array = {'jitter_PP_BOLD_PC_data', 'thres_PP_BOLD_PC_data'};

% PnP compairson
elseif isequal(jitter_data_type, 'PnP') && isequal(thres_data_type,'PnP')

    %Array of data types
    data_type_array = {'jitter_PnP_BOLD_PC_data', 'thres_PnP_BOLD_PC_data'};    
    
end

%Load voxelInfo for that ROI or masked
cd(ROI_dir)
load(fullfile(ROI_dir,'voxelInfo.mat'));

%Unique filenames for the constume built ROIs 
try
   
   %Define voxel coordinates
   vXYZ = voxelInfo;
   
 
end

try

    %Define voxel coordinates     
    vXYZ = xyz.vXYZ; 

end

%% Brain Subset Cut - Grid voxels defined by Find_adjacent_voxels.m

%If use grid voxels
if isequal(voxel_method, 2)
    
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
    
else
    
    disp('Using all voxels!')

end

%Specify image size (should be 91x109x91);
sizeMat = [91,109,91]; 

%Get the list of index of useful voxels sizeMat(2:4), is the range of X, Y, Z indices. sizeMat(1) is time
outInd = mat2ind(vXYZ,sizeMat);  %mat2ind creates an 1D index for the given matrix size (3D)

%Loop over data types (e.g., main data and subtraction data)
for type = 1:size(data_type_array, 2)

    %Define current data type
    current_data_type = data_type_array{type};
    
    disp(['Cutting voxels for ',current_data_type])
    
    %Initialize NaN matrix to store BOLD voxel values
    group_pc_voxel_subset_data = nan(size(outInd,1),size(eval(current_data_type),1),size(eval(current_data_type),3));

    %This for loop cycles through all of the participants
    for sess = 1:size(eval(current_data_type),3)

        %This for loop cycles through all time points per participant
        for timepoint = 1:size(eval(current_data_type),1) 

            %Reshape session voxel data to 91, 109, 91 for each time point and subject
            img = reshape(eval([current_data_type,'(timepoint,:,sess)']),[91, 109, 91]);
            
            %METHOD 1: use only single center voxel value
            if isequal(voxel_method, 1)

                %Select subset of voxels [voxels x time x subjects] - Grab the data from all the voxels specified in outInd form the 3D image
                group_pc_voxel_subset_data(:,timepoint,sess) = img(outInd);
            
            %METHOD 2: average signal from adjacent voxels
            elseif isequal(voxel_method, 2)
                
                %Loop over voxels
                for voxel = 1:size(vXYZ,2)
                    
                    %Find the mean BOLD signal for all voxels within spatial
                    %cluster, includes center and adjacent voxels - nanmean([adjacent voxels + center voxel])
                    group_pc_voxel_subset_data(voxel,timepoint,sess) = nanmean([img(all_adjacent_voxels_cell{voxel}); img(outInd(voxel))]);

                end
            
            end
            
        end

    end

    %Rename variable
    eval([current_data_type,' = group_pc_voxel_subset_data;']) % this is the ROI percent timecourse averaged across sessions
    
end

%% Subtraction and Baseline Datamean([img(all_adjacent_voxels_cell{voxel});img(outInd(voxel))])

disp('Baselining data from prestim period')

%Specify baseline period (defined between 1 and 41; stimulus onset at 21)
baseline_window = 1:20;

%Subtraction - PP minus PnP Trials
if isequal(jitter_data_type, 'PP_minus_PnP')
    
    %Rename
    PP_epochs_subjects_list = jitter_PP_epochs_subjects_list;
    PnP_epochs_subjects_list = jitter_PnP_epochs_subjects_list;
    group_PP_BOLD_PC_data = jitter_PP_BOLD_PC_data;
    group_PnP_BOLD_PC_data = jitter_PnP_BOLD_PC_data;
    
    %Remove subjects without PP and PnP trial instances (subject 603 does
    %not have PP trials)
    if length(PnP_epochs_subjects_list) > length(PP_epochs_subjects_list)
           
        group_PnP_BOLD_PC_data(:,:,not(ismember(PnP_epochs_subjects_list,PP_epochs_subjects_list))) = [];
        
        %Save list of the subjects tested
        tested_subjects_list = PnP_epochs_subjects_list(ismember(PnP_epochs_subjects_list,PP_epochs_subjects_list));
    
        %Check if there is PP subject but not in PnP (edge case)
        if size(group_PP_BOLD_PC_data,3) > size(group_PnP_BOLD_PC_data,3)
                       
           group_PP_BOLD_PC_data(:,:,not(ismember(PP_epochs_subjects_list,tested_subjects_list))) = [];
           
        end
        
    % If more PP than PnP subjects (Blank trials analysis many subjects
    % without PnP trials)
    elseif length(PP_epochs_subjects_list) > length(PnP_epochs_subjects_list)

        group_PP_BOLD_PC_data(:,:,not(ismember(PP_epochs_subjects_list,PnP_epochs_subjects_list))) = [];
        
        %Save list of the subjects tested
        tested_subjects_list = PP_epochs_subjects_list(ismember(PP_epochs_subjects_list,PnP_epochs_subjects_list));
        
    else
        
        %Create list of the subjects tested
        tested_subjects_list = PP_epochs_subjects_list;
    
    end
        
    %Subtract main data from subtract data
    jitter_group_pc_data = group_PP_BOLD_PC_data - group_PnP_BOLD_PC_data;
    
    %Rename
    jitter_tested_subject_list = tested_subjects_list;

end

%Subtraction - PP minus PnP Trials
if isequal(thres_data_type, 'PP_minus_PnP')
    
    %Rename
    PP_epochs_subjects_list = thres_PP_epochs_subjects_list;
    PnP_epochs_subjects_list = thres_PnP_epochs_subjects_list;
    group_PP_BOLD_PC_data = thres_PP_BOLD_PC_data;
    group_PnP_BOLD_PC_data = thres_PnP_BOLD_PC_data;
    
    %Remove subjects without PP and PnP trial instances (subject 603 does
    %not have PP trials)
    if length(PnP_epochs_subjects_list) > length(PP_epochs_subjects_list)
           
        group_PnP_BOLD_PC_data(:,:,not(ismember(PnP_epochs_subjects_list,PP_epochs_subjects_list))) = [];
        
        %Save list of the subjects tested
        tested_subjects_list = PnP_epochs_subjects_list(ismember(PnP_epochs_subjects_list,PP_epochs_subjects_list));
    
    % If more PP than PnP subjects (Blank trials analysis many subjects
    % without PnP trials)
    elseif length(PP_epochs_subjects_list) > length(PnP_epochs_subjects_list)

        group_PP_BOLD_PC_data(:,:,not(ismember(PP_epochs_subjects_list,PnP_epochs_subjects_list))) = [];
        
        %Save list of the subjects tested
        tested_subjects_list = PP_epochs_subjects_list(ismember(PP_epochs_subjects_list,PnP_epochs_subjects_list));
        
    else
        
        %Create list of the subjects tested
        tested_subjects_list = PP_epochs_subjects_list;
    
    end
        
    %Subtract main data from subtract data
    thres_group_pc_data = group_PP_BOLD_PC_data - group_PnP_BOLD_PC_data;
    
    %Rename
    thres_tested_subject_list = tested_subjects_list;

end

% %Jitter vs threshold
% %Remove subjects without PP and PnP trial instances
% if length(thres_tested_subject_list) > length(blank_tested_subject_list)

    thres_group_pc_data(:,:,not(ismember(thres_tested_subject_list, jitter_tested_subject_list))) = [];

    %Save list of the subjects tested
    tested_subjects_list_thres = thres_tested_subject_list(ismember(thres_tested_subject_list,jitter_tested_subject_list));
% 
% % If more PP than PnP subjects (Blank trials analysis many subjects
% % without PnP trials)
% elseif length(blank_tested_subject_list) > length(thres_tested_subject_list)

    jitter_group_pc_data(:,:,not(ismember(jitter_tested_subject_list,thres_tested_subject_list))) = [];

    %Save list of the subjects tested
    tested_subjects_list_jitter = jitter_tested_subject_list(ismember(jitter_tested_subject_list,thres_tested_subject_list));
% 
% end

%Check dimensions are correct
if ~isequal(size(thres_group_pc_data,3),size(jitter_group_pc_data,3))
    
    error('Data dimensions are mistmatched between jitter and threshold data')
    
elseif ~isequal(tested_subjects_list_jitter,tested_subjects_list_thres)
    
    error('Subject list mismatch between jitter and threshold data')
    
end

%Define direction of subtraction
if isequal(sub_direction,'No_Report')

    %Subtract main data from subtract data
    group_pc_data = thres_group_pc_data - jitter_group_pc_data;

elseif isequal(sub_direction,'ISI')
    
     %Subtract main data from subtract data
    group_pc_data = jitter_group_pc_data - thres_group_pc_data;
    
end

%Calculate baseline values [voxel x time x subjects]
group_pc_baseline = squeeze(nanmean(group_pc_data(:,baseline_window,:),2));

%Loop over time points
for time = 1:41

    %Subtract baseline from all voxels at all times
    group_pc_baselined_data(:,time,:) = squeeze(group_pc_data(:,time,:)) - group_pc_baseline;

end

%Clear unnecessary files
%clearvars group_pc_data group_pc_baseline group_PP_BOLD_PC_data group_PnP_BOLD_PC_data

%% Create Voxel Neigherhood - 2D Matrix

disp('Generating voxel proximity matrix')

%Setup empty hood matrix [total voxels x total voxels]
voxel_hood_matrix = zeros(size(vXYZ,2),size(vXYZ,2));

%Compare location (XYZ) for all voxels

%Loop over voxels
for test_voxel = 1:size(vXYZ,2)
    
        %Define test voxel coordinates (this voxel is compared to all others)
        test_xyz = vXYZ(:,test_voxel);
                    
    %Loop over voxels (again)
    for comp_voxel = 1:size(vXYZ,2)
        
        %Define compare voxel coordinates
        comp_xyz = vXYZ(:,comp_voxel);
        
        %Compare voxel dimension locations (XYZ) - each dimension must be within a specified range (distance threshold)
        if ~isempty(find(abs(test_xyz(1)-comp_xyz(1)) <= distance_threshold && abs(test_xyz(2)-comp_xyz(2)) <= distance_threshold && abs(test_xyz(3)-comp_xyz(3)) <= distance_threshold))

            %Fill matrix with adjacent voxel = 1
            voxel_hood_matrix(test_voxel,comp_voxel) = 1;
            voxel_hood_matrix(comp_voxel,test_voxel) = 1; %Add mirrored result
            
        end
    
    end
    
end

%Plot voxel hood (white = 1; black = 0)
%imshow(voxel_hood_matrix)

%% Run Permutation Test

disp('Running Permutation Test')

tic

%INPUTs: 3D data matrix, 2D Neighborhood matrix
[pval, t_orig, clust_info, seed_state, est_alpha, mn_clust_mass] = MR_clust_perm1_iceeg_sumt_rand(group_pc_baselined_data, voxel_hood_matrix, num_permutation, 0.05, 0, 0.05, 2, [], 0);

toc

%Save output
cd(fig_folder)
save('pc_cluster_sumt_rand_005.mat','pval','t_orig','clust_info','seed_state','est_alpha','mn_clust_mass' ...
,'group_pc_baselined_data', 'tested_subjects_list');

%% Plot brain map of cluster results images on MNI template

disp('Plotting brain maps')

%Load cluster stats
cd(fig_folder)
load('pc_cluster_sumt_rand_005.mat')

%Define the plotting method (1 = PC/BOLD; 2 = T-values; 3 = PC without stats threshold; 4 = Constant)
method_cell = {'PC_nostat','PC','t_values','constant'};

%Loop over methods
for plotting_method = 1:length(method_cell)
    
    %Define method names
    method_name = method_cell{plotting_method};
    
    disp(['Plotting ', method_name, ' maps'])

    %Figure folders
    fig_img_folder = fullfile(save_dir, ROI_name, ['Cluster_maps_',method_name]);
    hdr_img_folder = fullfile(save_dir, ROI_name, ['Source_imgs_',method_name]);

    %Create figure folders
    mkdir(fig_img_folder)
    mkdir(hdr_img_folder)

    %Enter figure folder
    cd(fig_folder)

    %Create empty matrix [Time x Voxels (91x109x91)]
    empty_matrix = nan(41,902629);

    %Average BOLD signal over subjects - Used in plotting percent change BOLD
    avg_group_pc_baselined_data = nanmean(group_pc_baselined_data,3);

    %Find significant clusters pvalue < 0.05 - Negative and positive clusters independently
    pos_clust = find(clust_info.pos_clust_pval < 0.05);
    neg_clust = find(clust_info.neg_clust_pval < 0.05);

    %If using only grid/subset voxels
    if isequal(voxel_method, 2)

        %METHOD: PC without significance thresholding
        if strcmp(method_name,'PC_nostat') %isequal(plotting_method, 3)

            %All grid voxels
            all_grid_voxels = 1:length(brain_grid_voxels);
           
            %Loop over time
            for time = 1:41
                
                %Find sig grid voxel t-values (sig voxels x time)
                sig_voxels_pc_values = avg_group_pc_baselined_data(all_grid_voxels,time);
            
                %Loop over significant voxels
                for voxel = 1:size(all_grid_voxels,2)

                    %Store PC in empty matrix
                    empty_matrix(time,outInd(all_grid_voxels(voxel))) = sig_voxels_pc_values(voxel);

                    %Store all adjacent voxels for each sig_grid_voxels in array
                    adj_voxels_array = cell2mat(all_adjacent_voxels_cell(voxel));

                     %Give adjacent voxels to post clusters constant value
                    empty_matrix(time,adj_voxels_array) = sig_voxels_pc_values(voxel);

                end   
            
           end
        
        % OTHER METHODS
        else
            
            %Loop over postive clusters
            for clust = 1:length(pos_clust)

                %Loop over time
                for time = 1:41

                    %Find postive significant cluster number/idx in voxel x time matrix
                    [sig_grid_voxels, col] = find(clust_info.pos_clust_ids(:,time) == pos_clust(clust));

                    %Include voxel adjacency in significant voxel designation (output = a cell of adjacent voxels for each of the significant voxels)
                    adj_voxels = all_adjacent_voxels_cell(sig_grid_voxels);

                    %METHOD: Percent change BOLD application for plotting
                    if strcmp(method_name,'PC') %isequal(plotting_method, 1)

                        %Find sig grid voxel t-values (sig voxels x time)
                        sig_voxels_pc_values = avg_group_pc_baselined_data(sig_grid_voxels,time);

                        %Loop over significant voxels
                        for voxel = 1:size(sig_grid_voxels,1)

                            %Store sig grid voxel t-values in empty matrix
                            empty_matrix(time,outInd(sig_grid_voxels(voxel))) = sig_voxels_pc_values(voxel);

                            %Store all adjacent voxels for each sig_grid_voxels in array
                            adj_voxels_array = cell2mat(adj_voxels(voxel));

                             %Give adjacent voxels to post clusters constant value
                            empty_matrix(time,adj_voxels_array) = sig_voxels_pc_values(voxel);

                        end

                    %METHOD: T-value application for plotting
                    elseif strcmp(method_name,'t_values') %isequal(plotting_method, 2)

                        %Find sig grid voxel t-values (sig voxels x time)
                        sig_voxels_t_values = t_orig(sig_grid_voxels,time);

                        %Loop over significant voxels
                        for voxel = 1:size(sig_grid_voxels,1)

                            %Store sig grid voxel t-values in empty matrix
                            empty_matrix(time,outInd(sig_grid_voxels(voxel))) = sig_voxels_t_values(voxel);

                            %Store all adjacent voxels for each sig_grid_voxels in array
                            adj_voxels_array = cell2mat(adj_voxels(voxel));

                             %Give adjacent voxels to post clusters constant value
                            empty_matrix(time,adj_voxels_array) = sig_voxels_t_values(voxel);

                        end

                    %METHOD: Constant value application for plotting
                    elseif strcmp(method_name,'constant') %isequal(plotting_method, 4)

                        %Initialize array variable used to covert adjacent cell to array
                        adj_voxels_array = [];

                        %Loop over sig_grid_voxels
                        for voxel = 1:size(sig_grid_voxels,1)

                            %Store all adjacent voxels for each sig_grid_voxels in array
                            adj_voxels_array = [adj_voxels_array; cell2mat(adj_voxels(voxel))];

                        end

                        %Apply constant for plotting.

                        %Give post clusters constant value - find the voxel all volume idx
                        empty_matrix(time,outInd(sig_grid_voxels)) = 0.5;

                        %Give adjacent voxels to post clusters constant value
                        empty_matrix(time,adj_voxels_array) = 0.5;

                    end

                end

            end

            %Loop over negative clusters
            for clust = 1:length(neg_clust)

                %Loop over time
                for time = 1:41

                    %Find negative significant cluster number/idx in voxel x time matrix
                    [sig_grid_voxels, col] = find(clust_info.neg_clust_ids(:,time) == neg_clust(clust));

                    %Include voxel adjacency in significant voxel designation (output = a cell of adjacent voxels for each of the significant voxels) 
                    adj_voxels = all_adjacent_voxels_cell(sig_grid_voxels);

                    %METHOD: Percent change BOLD application for plotting
                    if strcmp(method_name,'PC') %isequal(plotting_method, 1)

                        %Find sig grid voxel t-values (sig voxels x time)
                        sig_voxels_pc_values = avg_group_pc_baselined_data(sig_grid_voxels,time);

                        %Loop over significant voxels
                        for voxel = 1:size(sig_grid_voxels,1)

                            %Store sig grid voxel t-values in empty matrix
                            empty_matrix(time,outInd(sig_grid_voxels(voxel))) = sig_voxels_pc_values(voxel);

                            %Store all adjacent voxels for each sig_grid_voxels in array
                            adj_voxels_array = cell2mat(adj_voxels(voxel));

                             %Give adjacent voxels to negative clusters constant value
                            empty_matrix(time,adj_voxels_array) = sig_voxels_pc_values(voxel);

                        end

                    %METHOD: T-value application for plotting
                    elseif strcmp(method_name,'t_values') %isequal(plotting_method, 2)

                        %Find sig grid voxel t-values (sig voxels x time)
                        sig_voxels_t_values = t_orig(sig_grid_voxels,time);

                        %Loop over significant voxels
                        for voxel = 1:size(sig_grid_voxels,1)

                            %Store sig grid voxel t-values in empty matrix
                            empty_matrix(time,outInd(sig_grid_voxels(voxel))) = sig_voxels_t_values(voxel);

                            %Store all adjacent voxels for each sig_grid_voxels in array
                            adj_voxels_array = cell2mat(adj_voxels(voxel));

                             %Give adjacent voxels to post clusters constant value
                            empty_matrix(time,adj_voxels_array) = sig_voxels_t_values(voxel);

                        end 

                    %METHOD: Constant value application for plotting
                    elseif strcmp(method_name,'constant') %isequal(plotting_method, 4)

                        %Initialize array variable used to covert adjacent cell to array
                        adj_voxels_array = [];

                        %Loop over sig_grid_voxels
                        for voxel = 1:size(sig_grid_voxels,1)

                            %Store all adjacent voxels for each sig_grid_voxels in array
                            adj_voxels_array = [adj_voxels_array; cell2mat(adj_voxels(voxel))];

                        end

                        %Apply constant for plotting.

                        %Give negative clusters constant value - find the voxel all volume idx
                        empty_matrix(time,outInd(sig_grid_voxels)) = -0.5;

                        %Give adjacent voxels to negative clusters constant value
                        empty_matrix(time,adj_voxels_array) = -0.5;

                    end

                end

            end
                      
        end
                              
    end

    %Rename variable Dimensions: 41x902629
    mean_2D_PC_cut = empty_matrix;

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

    %Define plotting scale by method

    %Percent change with stats
    if strcmp(method_name,'PC') % isequal(plotting_method, 1) || isequal(plotting_method, 4)

        model.range(:,2) = [0.01; .2]; %Percent change threshold (Default = 0.15 to 1)
        model.range(:,3) = [0.01; .2]; %Percent change threshold (Default = 0.15 to 1)

    %Percent change without stats
    elseif strcmp(method_name,'PC_nostat')% isequal(plotting_method, 1) || isequal(plotting_method, 4)

        model.range(:,2) = [0.1; 1]; %Percent change threshold (Default = 0.15 to 1)
        model.range(:,3) = [0.1; 1]; %Percent change threshold (Default = 0.15 to 1)

    %T-values
    elseif strcmp(method_name,'t_values') %isequal(plotting_method, 2)

        model.range(:,2) = [0.01; 7];
        model.range(:,3) = [0.01; 7];

    %Constant scale
    elseif strcmp(method_name,'constant') %isequal(plotting_method, 3)

        model.range(:,2) = [0.15; 1]; 
        model.range(:,3) = [0.15; 1]; 

    end

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

end
