%% Percent Change Plot Script - UPDATED: Visual No Report v1.12

% Written by: Josh Ryu and Sharif I. Kronemer
% Date: 3/21/2017
% Edited: 12/4/2020

% Before running this script signal extraction should be completed
% This script completes the following functions in this order:
%(1) Applies a high pass filter 1/128 to extracted signal
%(2) Regresses out head motion 
%(3) Removes time points with FD and DVARS values above defined thresholds
%(4) Calculates a baseline value for extracted signal 
%(5) Calculate percent change

function Percent_Change_Map_Code_jitter_norm_new_bin_filter_Function(run_location, relevant_dir, ID, rootpath, run_num, run_list)

    tic
        
    %Loop over runs
    for run = 1:run_num

        %disp(['Running - Run ', num2str(run)])
        disp(['Running - ', run_list(run).name])

        %Subject specific run data folder
        %datafolder = fullfile(rootpath,['Run_', num2str(run)],'axialSlices'); 
        datafolder = fullfile(rootpath,run_list(run).name,'axialSlices'); 

        % Preprocessed image folder

        %Local
        if isequal(run_location, 'l')

            %rp_folder = fullfile('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Subject Analysis\',ID,'Perception Task',relevant_dir,'MRI Session','MRI Analysis',['Preprocessed Images\Run_',num2str(run)]);
            %rp_folder = fullfile('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Subject Analysis\',ID,'Perception Task',relevant_dir,'MRI Session','MRI Analysis','Preprocessed Images',run_list(run).name);
            rp_folder = fullfile('S:\HNCT No Report Paradigm\Subject Analysis MRI',ID,'Perception Task',relevant_dir,'MRI Session','MRI Analysis','Preprocessed Images',run_list(run).name);
            
        %Server
        elseif isequal(run_location, 's')

            %rp_folder = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Subject Analysis',ID,'Perception Task',relevant_dir,'MRI Session','MRI Analysis',['Preprocessed Images/Run_',num2str(run)]);
            %rp_folder = fullfile('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Subject Analysis',ID,'Perception Task',relevant_dir,'MRI Session','MRI Analysis','Preprocessed Images',run_list(run).name);
            rp_folder = fullfile('/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI',ID,'Perception Task',relevant_dir,'MRI Session','MRI Analysis','Preprocessed Images',run_list(run).name);

        end

        % Opens RP text file for motion regress step
        cd(rp_folder);
        rpfile_dir = dir('rp_*'); %Opens rp text file (generated during preprocessing)
        rp_data = readtable(rpfile_dir(1).name); % rp file name needs to be changed

        %% Apply a high pass filter to extracted voxel data

        % This implementation utilizes the frequency domain representation, and
        % deletes all the frequencies we are not interested in. 
        cd(datafolder)
        load('xyz_planes.mat'); 

        filtered_img = img;

        %Define image size depending on subject
        if isequal(ID, '458DW') && isequal(relevant_dir, 'Center Relevant') && isequal(run,1)

            %Specific case where wrong sequence was run for Run 1
            img_size = [600 91 109 91];

        else

            img_size = [700 91 109 91];

        end

        %Loop over dimensions
        for i = 1:img_size(2)% X dim
            for j = 1:img_size(3)% Y dim
                for kj = 1:img_size(4)% Z dim
                        Y = img(:,i,j,kj); % Y = a single voxel tracked through time
                    if isempty(nonzeros(Y))% For all zero's, skips to save time
                        filtered_img(:,i,j,kj) = Y;
                    else
                        k = length(Y);
                        HParam = 128; %frequency cut off
                        RT = 1;
                        n = fix(2*(k*RT)/HParam + 1); %Defines the number of regressions generated for filtering below 1/128
                        X0 = spm_dctmtx(k,n);%SPM function that makes a regression matrix to regress out frequencies below 1/128
                        S = X0(:,2:end); %Filtering matrix
                        IY = interp1(find(~isnan(Y)),Y(find(~isnan(Y))),[1:k]); %olate over nan values; makes a straight line between pt previous to nan and post nan.
                        IY = IY';
                        filtered_img(:,i,j,kj) = IY - S*(S'*IY); 
                    end
                end
            end
        end

        % Filtered image
        img = filtered_img; 

        %% Motion regressed out from extracted voxel data

        disp(['Regressing out head motion ....'])

        regressed_img = img;

        % Checks if motion data and img have the same dimension
        if size(rp_data,1) ~= size(img, 1)
            error('Motion data and image dimension mismatch')
        end

        % Loop over dimensions
        for i = 1:img_size(2)% X dim
            for j = 1:img_size(3)% Y dim
                for kj = 1:img_size(4)% Z dim
                    vox_tc = img(:,i,j,kj);
                    if isempty(nonzeros(vox_tc))%For all zero's, skips to save time
                        regressed_img(:,i,j,kj) = vox_tc;
                    else
                        mean_corrected = vox_tc - mean(vox_tc); %De-mean motion time course
                        [b, bint, r]= regress(mean_corrected, rp_data);%Uses matlab's regress function to regress out head motion
                        regressed_img(:,i,j,kj) = r + mean(vox_tc);%Re-mean motion time course
                    end
                end
            end
        end

        %% Outlier Rejection: DVARS and FD - Set Thresholds Here

        disp(['Applying outlier rejection criteria ....'])

        %DVARS
        % %To check if motion data and img have the same dimension
        % if length(rp_data) ~= size(img, 1)
        %     error('Motion data and image dimension mismatch')
        % end

        numVox = size(xyz.vXYZ,2);
        DVAR                    = DVARS(img,numVox);
        DVARthresh              = 5;%Enter DVAR threshold here (Default threshold = 5)
        DVAR_Ind                = find(DVAR>DVARthresh);

        %FID
        FIDthresh               = 0.3;%Enter FID threshold here (Default threshold = 0.3)
        [fd_Ind fDisp]          = moveArtCalc1(rp_data, FIDthresh);

        %Save variables
        cd(datafolder)
        save('artifact_timepoints_filtered.mat','DVARthresh', 'DVAR_Ind','FIDthresh', 'fd_Ind', 'regressed_img', '-v7.3')

        %% Calculate baseline of extracted voxel data

        % Baseline = the mean of the entire scan session (except for those times
        % points that did not pass the DVARS and FD thresholds above). 

        disp(['Calculating Baseline ....'])
        %The 'img' is the 4D image that is generated after extracting data for each voxel.
        %The baseline code below will find the average signal across time for each
        %voxel resulting in a 3D matrix with the dimentions 91 x 109 x 91.

        artifact_rm_img = regressed_img;
        artifact_rm_img(DVAR_Ind,:,:,:) = NaN; %At a particular time pt all voxels are rejected 
        %artifact_rm_img(DVAR_Ind+1,:,:,:) = NaN; %Added by josh
        artifact_rm_img(fd_Ind,:,:,:) = NaN; %At a particular time pt all voxels are rejected 
        %artifact_rm_img(fd_Ind+1,:,:,:) = NaN; %Added by josh
        %artifact_rm_img(roiSNR_Ind,:,:,:) = NaN; %Uncomment after implementing SNR

        extracted_signal_baseline = nanmean(permute(artifact_rm_img,[2 3 4 1]),4);%1st reorders the dimensions from [time X Y Z] -> [X Y Z time], then finds mean across the time dimension

        %Save variables
        cd(datafolder)
        %save('extracted_signal_baseline_filtered.mat','extracted_signal_baseline','artifact_rm_img','regressed_img', '-v7.3')
        save('extracted_signal_baseline_filtered.mat','extracted_signal_baseline', '-v7.3') %Only save baseline signal - not the img files

        %% Calculate Percent Change

        % Percent change equation = ((img data/img baseline)-1)*100

        disp(['Calculating Percent Change ....'])
        percent_change = nan(size(artifact_rm_img));%creates empty nan matrix

        for y = 1:size(artifact_rm_img,1)%Note that first dimension of this variable is time 

            artifact_rm_img_rs = reshape(artifact_rm_img(y,:,:,:),size(artifact_rm_img,2),size(artifact_rm_img,3),size(artifact_rm_img,4)); %transforms 4D to 3D matrix for a single time pt
            percent_change(y,:,:,:) = (artifact_rm_img_rs./extracted_signal_baseline -1)*100; %Percent change function

        end

        %Save variables
        cd(datafolder)
        save('extracted_signal_pc_filtered.mat','percent_change','-v7.3');

        clearvars img artifact_rm_img artifact_rm_img_rs extracted_signal_baseline
            
    end

    toc

end