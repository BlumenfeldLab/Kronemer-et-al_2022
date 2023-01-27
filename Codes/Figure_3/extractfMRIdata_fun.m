% batchfMRI_data: This function has three purposes:
% (1) Generate an spm design matrix for each run.  
%     This is treated as a dummy matrix that has no real purpose other than
%     to act as an input to marsbar functions.
% (2) Using marsbar, extract the ROI voxel values for every run. This
%     requires the use of an SPM design matrix, as mentioned in (1).
%     Typically, a whole brain gray matter mask is used as the ROI. 
% (3) From the marsbar output, format the fmri data outputs for each run
%     into a 4D matrix with the dimensions (time, x, y, z). 
% 
% INPUTS: 
%   resultsDir : directory to save results. 
%   aFile  : batch file that specifies which runs to use
%   dataDir: directory that contains patient fmri .swm files
%   roiPath: directory that contains mat file with ROI voxel info
%   isControl: 1 (for analyzing control data)

%Edited by: Sharif I. Kronemer
%Data edited: 12/4/2020
 
function extractfMRIdata_fun(resultsDir, dataDir, roiPath, run_location)

    tic

    % Revision History
    % originally main_estimation.m from xiaoxiao
    % change variable names such as matList, matDir since they reflect .img
    % files, not mat files
    % Convert to a function to specify different rois
    % ************** Set appropriate paths ******************

    % Marsbar path
    if isequal(run_location, 'l')

        a1 = genpath('Z:\Sharif Data F\HNCT fMRI Study\fMRI Analysis\Analysis Codes\marsbar-0.44');

    elseif isequal(run_location, 's')

        a1 = genpath('mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/marsbar-0.44');

    end

    %a2 = genpath('W:\HNCT fMRI Study\spm12_radiological');
    %a3 = genpath('/mnt/cae-lxs1/cae_code/fMRI libraries/spm8_mod');
    %b1 = genpath('/mnt/cae-lxs1/cae_code/CAE analysis code v1.8');
    %rmpath(a1) % Having marsbar on during GLM estimation will cause problems.
    %addpath(a2)
    %addpath(a3)
    %addpath(b1)

    % ******************************************************* 
    resultSPMPath   = fullfile(resultsDir,'SPM_data');
    resultMarsPath  = fullfile(resultsDir,'Mars_data');
    resultAxialPath = fullfile(resultsDir,'axialSlices');
    resultAxialPCPath = fullfile(resultsDir,'PCaxialSlices');
    
    % *************** Create SPM data ******************
    if ~isdir(resultSPMPath)
        mkdir(resultSPMPath)
    end
    
    cd(resultSPMPath)
    GLM_SPM12_estim(resultSPMPath,dataDir,1) %0.72 for HCP TR %1 for HNCTw
    
    % *************** Create Marsbar data ******************
    if ~isdir(resultMarsPath)
        mkdir(resultMarsPath)
    end
    
    addpath(a1) 

    %% need to change this. 
    spmFilePath = fullfile(resultSPMPath);
    %%

    auto_marsbar_spm12(spmFilePath,roiPath,resultMarsPath,1);

    % *************** Create Axial Slices ******************
    if ~isdir(resultAxialPath)
        mkdir(resultAxialPath)
        %rmdir(resultAxialPath, 's')
    end
    
    addpath(a1) 
    formatAxialSlices2(fullfile(resultMarsPath,'marsbar_data.mat'),resultAxialPath)

    %Save movement information from .ref file. 
    
    %% change this. 
    dataFolder = dataDir;
    
    moveFile = dir(fullfile(dataFolder,'rp_*'));
    moveFileName = moveFile.name;
    
    rp_data = load(fullfile(dataFolder,moveFileName));  
    save(fullfile(resultAxialPath,'xyz_planes.mat'),'rp_data','-append')

    % *******************************************************
    if exist(fullfile(resultMarsPath,'marsbar_data.mat'))~=2
        
        resultMarsPath
        
    end 

    toc

end
