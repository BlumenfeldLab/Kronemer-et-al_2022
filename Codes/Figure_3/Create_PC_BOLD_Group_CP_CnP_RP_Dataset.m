%% Create Group Percent Change BOLD Data Set - CP and CnP Trials Across Repor Paradigm

%The purpose of this code is to aggregate the percent change data from
%indivdual subjects and average to create a group percent change plot.

%Written by: Sharif I. Kronemer
%Date: 5/6/2021

clear

%% Prompts

%Run location
prompt1 = 'Run Location (l/s): ';
run_location = input(prompt1, 's');

%Post-stim trial duration
prompt2 = 'Post-stim duration (15/1): ';
poststim_delay = input(prompt2, 's');

%% Directories 

%Local
if isequal(run_location, 'l')
    
    %Add directories to path
    addpath(genpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\Extract voxel data\Supplementary functions'));
    addpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\spm12_radiological');
    addpath('S:\Data8\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\Behavioral Analysis')
    addpath('Z:\Sharif Data F\HNCT fMRI Study\fMRI Analysis\Analysis Codes\fMRI Behavioral Analysis')
    
    %Define template image
    template_image = 'S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\spm12_radiological\canonical\single_subj_T1.nii';
    
    %Group directory
    group_folder = 'Z:\Sharif Data F\HNCT fMRI Study\fMRI Analysis\Group Analysis\Percent Change Maps';
    
    %Subject direcotry 
    Report_subject_folder = 'Z:\Sharif Data F\HNCT fMRI Study\fMRI Analysis\Subject Analysis';
    
    %Save directory
    save_dir = 'Z:\Sharif Data F\HNCT fMRI Study\fMRI Analysis\Group Analysis\Group Data\BOLD PC';
    
%Server
elseif isequal(run_location, 's')
    
    %Add directories to path
    addpath(genpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Extract voxel data/Supplementary functions'));
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/spm12_radiological');
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/Behavioral Analysis')
    addpath('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/fMRI Behavioral Analysis')
    
    %Define template image
    template_image = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/spm12_radiological/canonical/single_subj_T1.nii';
    
    %Group directory
    group_folder = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Group Analysis/Percent Change Maps'; 
    
    %Subject directory
    Report_subject_folder = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Subject Analysis';
    
    %Save directory
    save_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Group Data/BOLD PC';
    
end

%% Group Variables and Subject Lists

%Find all subject folders
Report_subject_list = dir(Report_subject_folder);
Report_subject_list = {Report_subject_list.name}';

%Setup subject count
subject_count = 0;

%Initialized variables
group_CP_BOLD_PC_data = [];
group_CnP_BOLD_PC_data = [];
all_included_subjects_list = {};

%% REPORT ONLY PARADIGM SUBJECTS - Aggregate percent change over subjects

%Loop over subject
for sub = 1:length(Report_subject_list) 
    
    %Select ID
    ID = Report_subject_list{sub};
           
    disp(['**Added Report Subject ',num2str(ID),'**'])

    %Rejecting subjects by behavioral performance
    bad_subject_idx = report_subject_rejection_by_behavior_MRI(ID);

    %Check if bad subject by behavior
    if isequal(bad_subject_idx, 1)

        disp(['Bad behavior - Skipping ', num2str(ID)])
        continue

    end
    
    %Special case subject rejection - 252NT Run 1,2,3 are rejected by
    %motion artifact and the remaining Run 4 has no CP 15s post stim
    %trials, so all permutations statistics are compromised
    if isequal(ID, '252NT')
        
        disp(['Bad motion - Skipping ', num2str(ID)])
        continue
       
    end 

    %Subject percent change data path 
    if isequal(run_location, 'l')

        rootpath = fullfile('Z:\Sharif Data F\HNCT fMRI Study\fMRI Analysis\Subject Analysis',ID,'MRI Analysis\Extracted voxel data\session_percent_change');

    elseif isequal(run_location, 's')

        rootpath = fullfile('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Subject Analysis',ID,'MRI Analysis/Extracted voxel data/session_percent_change');

    end

    %Check if root path exists
    if exist(rootpath)

        %Check if the mat file exists
        if ~exist(fullfile(rootpath,'extracted_signal_PC_CP15_norm_binned_cuts_filtered.mat')) ...
           && ~exist(fullfile(rootpath,'extracted_signal_PC_CnP15_norm_binned_cuts_filtered.mat'))

            disp(['Cannot find data - Skipping ', num2str(ID)])
            continue

        end

    else

       disp(['Cannot find directory - Skipping ', num2str(ID)])
       continue

    end
            
    %Add subject ID to subject list
    all_included_subjects_list = [all_included_subjects_list; ID];
       
    %Load CP percent change data
    cd(rootpath)
    load(['extracted_signal_PC_CP',num2str(poststim_delay),'_norm_binned_cuts_filtered.mat'])

    % Create group matrix of time points x voxels x subjects
    if isempty(group_CP_BOLD_PC_data) 

        group_CP_BOLD_PC_data(:,:,1)  = subject_img_3D_PC_cut(:,:,:); 

    else

        group_CP_BOLD_PC_data(:,:,(end+1))  = subject_img_3D_PC_cut(:,:,:);

    end 

    %Load CnP percent change data
    cd(rootpath)
    
    %Some CnP variables have a slightly different naming scheme
    try
        load(['extracted_signal_PC_CnP',num2str(poststim_delay),'_norm_binned_cuts_filtered.mat'])

    catch
        load(['extracted_signal_PC_CnP',num2str(poststim_delay),'s_norm_binned_cuts_filtered.mat'])

    end
    
    % Create group matrix of time points x voxels x subjects
    if isempty(group_CnP_BOLD_PC_data) 

        group_CnP_BOLD_PC_data(:,:,1)  = subject_img_3D_PC_cut(:,:,:); 

    else

        group_CnP_BOLD_PC_data(:,:,(end+1))  = subject_img_3D_PC_cut(:,:,:);

    end

end

%% Save Report and No Report Paradigm Dataset Matrices

%Save subject BOLD matrices
cd(save_dir)
save(['Group_PC_BOLD_CP_CnP_poststim_',num2str(poststim_delay),'s_data.mat'], 'group_CP_BOLD_PC_data', 'group_CnP_BOLD_PC_data', 'all_included_subjects_list', '-v7.3')

%{
%% Plot percent change images on MNI template

%Name trial types
trial_types_cell = {'CP','CnP','CP_minus_CnP'};

%Loop over trial types
for type = 1:length(trial_types_cell)

    %Define current trial type
    trial_type = trial_types_cell{type};
    
    disp(['Plotting BOLD PC Maps - ',trial_type])
    
    %% Average BOLD across subjects

    %CP Trials
    if isequal(trial_type, 'CP')

        %Average Across Subjects
        mean_PC_data = nanmean(group_CP_BOLD_PC_data,3);

    %CnP Trials
    elseif isequal(trial_type, 'CnP')

        %Average Across Subjects
        mean_PC_data = nanmean(group_CnP_BOLD_PC_data,3);

    %CP-CnP Data
    elseif isequal(trial_type, 'CP_minus_CnP')

        %Check if PC data variables are the same size
        if ~isequal(size(group_CP_BOLD_PC_data,3),size(group_CnP_BOLD_PC_data,3))

            error('Different sized PC data to substract!') 

        end

        %Subtract mean CP and CnP resposnes
        mean_PC_data = nanmean(group_CP_BOLD_PC_data,3) - nanmean(group_CnP_BOLD_PC_data,3);
        %mean_PC_data1 = nanmean((group_CP_BOLD_PC_data-group_CnP_BOLD_PC_data),3);

    end    

    %Define final image folders
    final_data_folder = fullfile(group_folder,['PC_maps_group_avg_',trial_type]);
    final_img_folder = fullfile(final_data_folder,'PC_act_deact_maps');

    %Make percent change map directories    
    mkdir(final_data_folder);
    mkdir(final_img_folder);

    %Begin plotting 
    cd(final_data_folder)

    %Loop over time points
    for time = 1:size(mean_PC_data,1)

        %Reshape into 91x109x91 matrix
        img_3D_pertime = reshape(mean_PC_data(time,:),[91, 109, 91]);

        data_pos{time} = zeros(size(img_3D_pertime));
        data_neg{time} = zeros(size(img_3D_pertime));
        data_pos{time} = (img_3D_pertime > 0.01).*img_3D_pertime; %Percent change values above 0.01 are consider positive
        data_neg{time} = (img_3D_pertime < -0.01).*img_3D_pertime; %Percent cahnge values below -0.01 are considered negative

    end

    %Open template voxelInfo
    if isequal(run_location, 'l')

        load('Z:\Sharif Data F\HNCT fMRI Study\fMRI Analysis\Subject Analysis\191AM\MRI Analysis\Extracted voxel data\Run 1 Movie\axialSlices\voxelInfo.mat')

    elseif isequal(run_location, 's')

        load('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Subject Analysis/191AM/MRI Analysis/Extracted voxel data/Run 1 Movie/axialSlices/voxelInfo.mat')

    end

    xyz = xyz.vXYZ;
    avw = load_nii(template_image);

    for m = 1:size(data_pos, 2)

        newimg = double(data_pos{m});
        avw.img = avw.img*0;

        for i = 1: size(xyz, 2)
            avw.img(xyz(1,i), xyz(2, i), xyz(3, i)) = newimg(size(newimg,1)-xyz(1,i)+1,xyz(2, i),xyz(3, i)); %flip L-R
        end

        avw.hdr.dime.datatype = 16;
        avw.hdr.dime.bitpix = int16(32);
        avw.hdr.dime.cal_max = max(max(max(avw.img)));
        avw.hdr.dime.roi_scale = 1;
        save_nii(avw, ['Act_time_' num2str(m)]);
        newimg = double(data_neg{m});
        avw.img = avw.img*0;

        for i = 1: size(xyz, 2)
            avw.img(xyz(1,i), xyz(2, i), xyz(3, i)) = -newimg(size(newimg,3)-xyz(1,i)+1, xyz(2, i), xyz(3, i)); %flip L-R
        end

        avw.hdr.dime.datatype = 16;
        avw.hdr.dime.bitpix = int16(32);
        avw.hdr.dime.cal_max = max(max(max(avw.img)));
        avw.hdr.dime.roi_scale = 1;
        save_nii(avw, ['Dec_time_' num2str(m)]);

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
    model.range(:,2) = [0.1; 0.5]; %Percent change threshold (Default = 0.15 to 1)
    model.range(:,3) = [0.1; 0.5]; %Percent change threshold (Default = 0.15 to 1)
    model.transform = 'axial'; %'axial','coronal','sagittal'
    model.axialslice = [-56:8:86]; %Defines the # of slices and slice thickness plotted; full range:[-72:2:102]
    model.coronalslice = [-92:8:52]; %Defines the # of slices and slice thickness plotted; full range:[-72:2:102]

    %Loop over time points
    for m = 1:size(data_pos, 2)

        %Setup model
        model.imgs{1, 1} = template_image;
        model.imgs{2, 1} = fullfile(final_data_folder, ['Act_time_' num2str(m) '.img']);
        model.imgs{3, 1} = fullfile(final_data_folder, ['Dec_time_' num2str(m) '.img']);

        %Run image generation
        display_slices_bai_1;

        cd(final_img_folder);
        set(gcf,'PaperPositionMode','auto');
        print('-dtiff', ['Time_' num2str(m)]);
        cd('..');
        close all;

    end 

end
%}