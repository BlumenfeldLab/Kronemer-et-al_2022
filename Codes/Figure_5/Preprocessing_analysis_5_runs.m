%% MRI Preprocessing - Raw Nifiti Files to Smoothed Nifiti

%The purpose of the code is to read raw nifiti files and run preprocessing.
%You need to reorient the MPRAGE image before running this code. Before
%running this code make sure to have reoriented the MPRAGE for each study
%session to the anterior commissure. 

% (1) Realign estimate
% (2) Realign reslice
% (3) Coregistration
% (4) Normalization
% (5) Smooth

%Written by: Sharif I. Kronemer
%Date: 3/21/2020

clear

%% Define parameters

% Select Subject Data to Analyze
prompt = 'Subject ID: '; %prompt ID
ID = input(prompt, 's');

% Location Relevant
prompt_2 = 'Relevant location (q/c): '; %prompt ID
relevant = input(prompt_2, 's');

% Define relevant directory
if isequal(relevant, 'q')
    
   relevant_location = 'Quadrant Relevant';
   
elseif isequal(relevant, 'c')
    
   relevant_location = 'Center Relevant';

end

% MR Number
prompt_3 = 'MRI number: '; %prompt ID
MRI_num = input(prompt_3, 's');

% Run code location
prompt_4 = 'Run code location (l/s): '; %prompt ID
run_location = input(prompt_4, 's');

%% Define Directories

%Local
if isequal(run_location, 'l')%r)
    
    %Add spm to path
    addpath('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\spm12_radiological')

    %Reoriented structural image directory
    %reoriented_structural = dir(['S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Subject Analysis\', ID, '\Perception Task\', relevant_location, '\MRI Session\MRI Analysis\spa*.nii']);
    reoriented_structural = dir(['Y:\HNCT No Report Paradigm\Subject Analysis MRI\', ID, '\Perception Task\', relevant_location, '\MRI Session\MRI Analysis\spa*.nii']);
   
    %Create exception for pb MR ID file names
    if isequal(size(reoriented_structural,1), 0)
        
        %reoriented_structural = dir(['S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Subject Analysis\', ID, '\Perception Task\', relevant_location, '\MRI Session\MRI Analysis\spb*.nii']);       
        reoriented_structural = dir(['Y:\HNCT No Report Paradigm\Subject Analysis MRI\', ID, '\Perception Task\', relevant_location, '\MRI Session\MRI Analysis\spb*.nii']);       

    end
    
    %Create reoriented image directory
    reoriented_structural = fullfile(reoriented_structural.folder, reoriented_structural.name);

    %Nii directory 
    %data_directory = ['S:\HNCT_NRP_Study\HNCT No Report Paradigm\Subject Raw Data\', ID, '\Perception Task\',relevant_location,'\MRI Session\MRI Data\Converted nii'];
    data_directory = ['M:\Subject Raw Data\', ID, '\Perception Task\',relevant_location,'\MRI Session\MRI Data\Converted nii'];

    %Save directory
    %save_directory = ['S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Subject Analysis\',ID,'\Perception Task\',relevant_location,'\MRI Session\MRI Analysis\Preprocessed Images'];
    save_directory = ['Y:\HNCT No Report Paradigm\Subject Analysis MRI\',ID,'\Perception Task\',relevant_location,'\MRI Session\MRI Analysis\Preprocessed Images'];

    %MNI directory 
    MNI_directory = 'S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Analysis Code\MRI Analysis\spm12_radiological\tpm\TPM.nii';

%Server
elseif isequal(run_location, 's')

    %Add spm to path
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/spm12_radiological')

    %Reoriented structural image directory
    %reoriented_structural = dir(['/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Subject Analysis/', ID, '/Perception Task/', relevant_location, '/MRI Session/MRI Analysis/s',MRI_num,'*.nii']);
    %reoriented_structural = dir(['/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Subject Analysis/', ID, '/Perception Task/', relevant_location, '/MRI Session/MRI Analysis/s*.nii']);
    reoriented_structural = dir(['/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI/', ID, '/Perception Task/', relevant_location, '/MRI Session/MRI Analysis/s*.nii']);
   
    %Create exception for pb MR ID file names
    if isequal(size(reoriented_structural,1), 0)
        
        %reoriented_structural = dir(['/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Subject Analysis/', ID, '/Perception Task/', relevant_location, '/MRI Session/MRI Analysis/spb*.nii']);       
        reoriented_structural = dir(['/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI/', ID, '/Perception Task/', relevant_location, '/MRI Session/MRI Analysis/spb*.nii']);       
    
    end
    
    %Create reoriented image directory
    reoriented_structural = fullfile(reoriented_structural.folder, reoriented_structural.name);
    
    %Nii directory 
    %data_directory = ['/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Subject Raw Data/', ID, '/Perception Task/',relevant_location,'/MRI Session/MRI Data/Converted nii'];
    data_directory = ['/mnt/Data15/Subject Raw Data/', ID, '/Perception Task/',relevant_location,'/MRI Session/MRI Data/Converted nii'];

    %Save directory
    %save_directory = ['/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Subject Analysis/',ID,'/Perception Task/',relevant_location,'/MRI Session/MRI Analysis/Preprocessed Images'];
    save_directory = ['/mnt/Data18/HNCT No Report Paradigm/Subject Analysis MRI/',ID,'/Perception Task/',relevant_location,'/MRI Session/MRI Analysis/Preprocessed Images'];
    
    %MNI directory 
    MNI_directory = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/spm12_radiological/tpm/TPM.nii';

end

%% Open SPM

spm fmri

%% Select nii files to preprocess

%Find all nii files

%NOTE: Commented out sections below manage unique subject session numbers. The default 
%values for 5 runs of the no-report paradgim are session numbers 4-8. However, when another 
%sequence was added (e.g., restarting a sequence in the study), the sequence numbers can change 
%and this section of the script must be updated so that the correct sequences are found. 

if isequal(ID, '571') && isequal(relevant, 'c')

    %Note: Session 7 is neglected because the participant hit the panic
    %button and run was restarted
    Run_1_nii = dir(fullfile(data_directory,['f',MRI_num,'-0004*']));
    Run_2_nii = dir(fullfile(data_directory,['f',MRI_num,'-0005*']));
    Run_3_nii = dir(fullfile(data_directory,['f',MRI_num,'-0006*']));
    Run_4_nii = dir(fullfile(data_directory,['f',MRI_num,'-0008*']));
    Run_5_nii = dir(fullfile(data_directory,['f',MRI_num,'-0009*']));
 
else

    %If the first run session number is 5
    Run_1_nii = dir(fullfile(data_directory,['f',MRI_num,'-0005*']));
    Run_2_nii = dir(fullfile(data_directory,['f',MRI_num,'-0006*']));
    Run_3_nii = dir(fullfile(data_directory,['f',MRI_num,'-0007*']));
    Run_4_nii = dir(fullfile(data_directory,['f',MRI_num,'-0008*']));
    Run_5_nii = dir(fullfile(data_directory,['f',MRI_num,'-0009*']));
     
%     %If the first run session number is 6
%     Run_1_nii = dir(fullfile(data_directory,['f',MRI_num,'-0006*']));
%     Run_2_nii = dir(fullfile(data_directory,['f',MRI_num,'-0007*']));
%     Run_3_nii = dir(fullfile(data_directory,['f',MRI_num,'-0008*']));
%     Run_4_nii = dir(fullfile(data_directory,['f',MRI_num,'-0009*']));
%     Run_5_nii = dir(fullfile(data_directory,['f',MRI_num,'-0010*']));
    
%     %If the first run session number is 4
%     Run_1_nii = dir(fullfile(data_directory,['f',MRI_num,'-0004*']));
%     Run_2_nii = dir(fullfile(data_directory,['f',MRI_num,'-0005*']));
%     Run_3_nii = dir(fullfile(data_directory,['f',MRI_num,'-0006*']));
%     Run_4_nii = dir(fullfile(data_directory,['f',MRI_num,'-0007*']));
%     Run_5_nii = dir(fullfile(data_directory,['f',MRI_num,'-0008*']));

%     %If the first run session number is 4,5,7,8,9
%     Run_1_nii = dir(fullfile(data_directory,['f',MRI_num,'-0004*']));
%     Run_2_nii = dir(fullfile(data_directory,['f',MRI_num,'-0005*']));
%     Run_3_nii = dir(fullfile(data_directory,['f',MRI_num,'-0007*']));
%     Run_4_nii = dir(fullfile(data_directory,['f',MRI_num,'-0008*']));
%     Run_5_nii = dir(fullfile(data_directory,['f',MRI_num,'-0009*']));

end

%Find file path of all nii files

%Run 1

%Initialize variable
Run_1_nii_filepath = {};

%Loop over rows
for row = 1:size(Run_1_nii, 1)
    
    %Store paths 
    Run_1_nii_filepath{row, 1} = fullfile(Run_1_nii(row).folder, Run_1_nii(row).name);
    
end

%Run 2

%Initialize variable
Run_2_nii_filepath = {};

%Loop over rows
for row = 1:size(Run_2_nii, 1)
    
    %Store paths 
    Run_2_nii_filepath{row, 1} = fullfile(Run_2_nii(row).folder, Run_2_nii(row).name);
    
end

%Run 3

%Initialize variable
Run_3_nii_filepath = {};

%Loop over rows
for row = 1:size(Run_3_nii, 1)
    
    %Store paths 
    Run_3_nii_filepath{row, 1} = fullfile(Run_3_nii(row).folder, Run_3_nii(row).name);
    
end

%Run 4

%Initialize variable
Run_4_nii_filepath = {};

%Loop over rows
for row = 1:size(Run_4_nii, 1)
    
    %Store paths 
    Run_4_nii_filepath{row, 1} = fullfile(Run_4_nii(row).folder, Run_4_nii(row).name);
    
end

%Run 5 

%Initialize variable
Run_5_nii_filepath = {};

%Loop over rows
for row = 1:size(Run_5_nii, 1)
    
    %Store paths 
    Run_5_nii_filepath{row, 1} = fullfile(Run_5_nii(row).folder, Run_5_nii(row).name);
    
end

tic

%% Realign and estimate

%Select the files of interest
matlabbatch{1}.spm.spatial.realign.estimate.data = {Run_1_nii_filepath, Run_2_nii_filepath, Run_3_nii_filepath, Run_4_nii_filepath, Run_5_nii_filepath};

matlabbatch{1}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estimate.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estimate.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estimate.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realign.estimate.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estimate.eoptions.weight = '';

%% Realign and write 

% Session 1
matlabbatch{2}.spm.spatial.realign.write.data(1) = cfg_dep('Realign: Estimate: Realigned Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','cfiles'));
matlabbatch{2}.spm.spatial.realign.write.roptions.which = [2 1];
matlabbatch{2}.spm.spatial.realign.write.roptions.interp = 4;
matlabbatch{2}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.realign.write.roptions.mask = 1;
matlabbatch{2}.spm.spatial.realign.write.roptions.prefix = 'r';

% Session 2
matlabbatch{3}.spm.spatial.realign.write.data(1) = cfg_dep('Realign: Estimate: Realigned Images (Sess 2)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','cfiles'));
matlabbatch{3}.spm.spatial.realign.write.roptions.which = [2 1];
matlabbatch{3}.spm.spatial.realign.write.roptions.interp = 4;
matlabbatch{3}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
matlabbatch{3}.spm.spatial.realign.write.roptions.mask = 1;
matlabbatch{3}.spm.spatial.realign.write.roptions.prefix = 'r';

% Session 3
matlabbatch{4}.spm.spatial.realign.write.data(1) = cfg_dep('Realign: Estimate: Realigned Images (Sess 3)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{3}, '.','cfiles'));
matlabbatch{4}.spm.spatial.realign.write.roptions.which = [2 1];
matlabbatch{4}.spm.spatial.realign.write.roptions.interp = 4;
matlabbatch{4}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
matlabbatch{4}.spm.spatial.realign.write.roptions.mask = 1;
matlabbatch{4}.spm.spatial.realign.write.roptions.prefix = 'r';

% Session 4
matlabbatch{5}.spm.spatial.realign.write.data(1) = cfg_dep('Realign: Estimate: Realigned Images (Sess 4)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{4}, '.','cfiles'));
matlabbatch{5}.spm.spatial.realign.write.roptions.which = [2 1];
matlabbatch{5}.spm.spatial.realign.write.roptions.interp = 4;
matlabbatch{5}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
matlabbatch{5}.spm.spatial.realign.write.roptions.mask = 1;
matlabbatch{5}.spm.spatial.realign.write.roptions.prefix = 'r';

% Session 5
matlabbatch{6}.spm.spatial.realign.write.data(1) = cfg_dep('Realign: Estimate: Realigned Images (Sess 5)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{5}, '.','cfiles'));
matlabbatch{6}.spm.spatial.realign.write.roptions.which = [2 1];
matlabbatch{6}.spm.spatial.realign.write.roptions.interp = 4;
matlabbatch{6}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
matlabbatch{6}.spm.spatial.realign.write.roptions.mask = 1;
matlabbatch{6}.spm.spatial.realign.write.roptions.prefix = 'r';

%% Coregistration

% Session 1
matlabbatch{7}.spm.spatial.coreg.estwrite.ref = {reoriented_structural};
matlabbatch{7}.spm.spatial.coreg.estwrite.source(1) = cfg_dep('Realign: Reslice: Mean Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
matlabbatch{7}.spm.spatial.coreg.estwrite.other(1) = cfg_dep('Realign: Reslice: Resliced Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));
matlabbatch{7}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{7}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{7}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{7}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{7}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{7}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{7}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{7}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

% Session 2
matlabbatch{8}.spm.spatial.coreg.estwrite.ref = {reoriented_structural};
matlabbatch{8}.spm.spatial.coreg.estwrite.source(1) = cfg_dep('Realign: Reslice: Mean Image', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
matlabbatch{8}.spm.spatial.coreg.estwrite.other(1) = cfg_dep('Realign: Reslice: Resliced Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));
matlabbatch{8}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{8}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{8}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{8}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{8}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{8}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{8}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{8}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

% Session 3
matlabbatch{9}.spm.spatial.coreg.estwrite.ref = {reoriented_structural};
matlabbatch{9}.spm.spatial.coreg.estwrite.source(1) = cfg_dep('Realign: Reslice: Mean Image', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
matlabbatch{9}.spm.spatial.coreg.estwrite.other(1) = cfg_dep('Realign: Reslice: Resliced Images', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));
matlabbatch{9}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{9}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{9}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{9}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{9}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{9}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{9}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{9}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

% Session 4
matlabbatch{10}.spm.spatial.coreg.estwrite.ref = {reoriented_structural};
matlabbatch{10}.spm.spatial.coreg.estwrite.source(1) = cfg_dep('Realign: Reslice: Mean Image', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
matlabbatch{10}.spm.spatial.coreg.estwrite.other(1) = cfg_dep('Realign: Reslice: Resliced Images', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));
matlabbatch{10}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{10}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{10}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{10}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{10}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{10}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{10}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{10}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

% Session 5
matlabbatch{11}.spm.spatial.coreg.estwrite.ref = {reoriented_structural};
matlabbatch{11}.spm.spatial.coreg.estwrite.source(1) = cfg_dep('Realign: Reslice: Mean Image', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
matlabbatch{11}.spm.spatial.coreg.estwrite.other(1) = cfg_dep('Realign: Reslice: Resliced Images', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));
matlabbatch{11}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{11}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{11}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{11}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{11}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{11}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{11}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{11}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

%% Normalization

% Session 1
matlabbatch{12}.spm.spatial.normalise.estwrite.subj(1).vol = {reoriented_structural};
matlabbatch{12}.spm.spatial.normalise.estwrite.subj(1).resample(1) = cfg_dep('Coregister: Estimate & Reslice: Resliced Images', substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));

% Session 2
matlabbatch{12}.spm.spatial.normalise.estwrite.subj(2).vol = {reoriented_structural};
matlabbatch{12}.spm.spatial.normalise.estwrite.subj(2).resample(1) = cfg_dep('Coregister: Estimate & Reslice: Resliced Images', substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));

% Session 3
matlabbatch{12}.spm.spatial.normalise.estwrite.subj(3).vol = {reoriented_structural};
matlabbatch{12}.spm.spatial.normalise.estwrite.subj(3).resample(1) = cfg_dep('Coregister: Estimate & Reslice: Resliced Images', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));

% Session 4
matlabbatch{12}.spm.spatial.normalise.estwrite.subj(4).vol = {reoriented_structural};
matlabbatch{12}.spm.spatial.normalise.estwrite.subj(4).resample(1) = cfg_dep('Coregister: Estimate & Reslice: Resliced Images', substruct('.','val', '{}',{10}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));

% Session 5
matlabbatch{12}.spm.spatial.normalise.estwrite.subj(5).vol = {reoriented_structural};
matlabbatch{12}.spm.spatial.normalise.estwrite.subj(5).resample(1) = cfg_dep('Coregister: Estimate & Reslice: Resliced Images', substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));

matlabbatch{12}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{12}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{12}.spm.spatial.normalise.estwrite.eoptions.tpm = {MNI_directory};
matlabbatch{12}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
matlabbatch{12}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{12}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
matlabbatch{12}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
matlabbatch{12}.spm.spatial.normalise.estwrite.woptions.bb = [-90 -126 -72
                                                              90 90 108];
matlabbatch{12}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
matlabbatch{12}.spm.spatial.normalise.estwrite.woptions.interp = 4;
matlabbatch{12}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';

%% Smooth

% Session 1
matlabbatch{13}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Estimate & Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{12}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{13}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{13}.spm.spatial.smooth.dtype = 0;
matlabbatch{13}.spm.spatial.smooth.im = 0;
matlabbatch{13}.spm.spatial.smooth.prefix = 's';

% Session 2
matlabbatch{14}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Estimate & Write: Normalised Images (Subj 2)', substruct('.','val', '{}',{12}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{2}, '.','files'));
matlabbatch{14}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{14}.spm.spatial.smooth.dtype = 0;
matlabbatch{14}.spm.spatial.smooth.im = 0;
matlabbatch{14}.spm.spatial.smooth.prefix = 's';

% Session 3
matlabbatch{15}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Estimate & Write: Normalised Images (Subj 3)', substruct('.','val', '{}',{12}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{3}, '.','files'));
matlabbatch{15}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{15}.spm.spatial.smooth.dtype = 0;
matlabbatch{15}.spm.spatial.smooth.im = 0;
matlabbatch{15}.spm.spatial.smooth.prefix = 's';

% Session 4
matlabbatch{16}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Estimate & Write: Normalised Images (Subj 4)', substruct('.','val', '{}',{12}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{4}, '.','files'));
matlabbatch{16}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{16}.spm.spatial.smooth.dtype = 0;
matlabbatch{16}.spm.spatial.smooth.im = 0;
matlabbatch{16}.spm.spatial.smooth.prefix = 's';

% Session 5
matlabbatch{17}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Estimate & Write: Normalised Images (Subj 5)', substruct('.','val', '{}',{12}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{5}, '.','files'));
matlabbatch{17}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{17}.spm.spatial.smooth.dtype = 0;
matlabbatch{17}.spm.spatial.smooth.im = 0;

%% RUN

spm_jobman('run', matlabbatch)

%% Move and delete files

disp('Moving preprocessed files')

%Make directories
mkdir(fullfile(save_directory,'Run_1'))
mkdir(fullfile(save_directory,'Run_2'))
mkdir(fullfile(save_directory,'Run_3'))
mkdir(fullfile(save_directory,'Run_4'))
mkdir(fullfile(save_directory,'Run_5'))

cd(data_directory)

%NOTE: Commented out sections below manage unique subject session numbers. The default 
%values for 5 runs of the no-report paradgim are session numbers 4-8. However, when another 
%sequence was added (e.g., restarting a sequence in the study), the sequence numbers can change 
%and this section of the script must be updated so that the correct sequences are found. 

%{
%Move smoothed files
movefile(['swrrf',MRI_num,'-0004*'], fullfile(save_directory,'Run_1'))
movefile(['swrrf',MRI_num,'-0005*'], fullfile(save_directory,'Run_2'))
movefile(['swrrf',MRI_num,'-0006*'], fullfile(save_directory,'Run_3'))
movefile(['swrrf',MRI_num,'-0007*'], fullfile(save_directory,'Run_4'))
movefile(['swrrf',MRI_num,'-0008*'], fullfile(save_directory,'Run_5'))

%Move rp text files
movefile(['rp_f',MRI_num,'-0004*'], fullfile(save_directory,'Run_1'))
movefile(['rp_f',MRI_num,'-0005*'], fullfile(save_directory,'Run_2'))
movefile(['rp_f',MRI_num,'-0006*'], fullfile(save_directory,'Run_3'))
movefile(['rp_f',MRI_num,'-0007*'], fullfile(save_directory,'Run_4'))
movefile(['rp_f',MRI_num,'-0008*'], fullfile(save_directory,'Run_5'))

%Move mean nii files
movefile(['meanf',MRI_num,'-0004*'], fullfile(save_directory,'Run_1'))
movefile(['meanf',MRI_num,'-0005*'], fullfile(save_directory,'Run_2'))
movefile(['meanf',MRI_num,'-0006*'], fullfile(save_directory,'Run_3'))
movefile(['meanf',MRI_num,'-0007*'], fullfile(save_directory,'Run_4'))
movefile(['meanf',MRI_num,'-0008*'], fullfile(save_directory,'Run_5'))

movefile(['rmeanf',MRI_num,'-0004*'], fullfile(save_directory,'Run_1'))
movefile(['rmeanf',MRI_num,'-0005*'], fullfile(save_directory,'Run_2'))
movefile(['rmeanf',MRI_num,'-0006*'], fullfile(save_directory,'Run_3'))
movefile(['rmeanf',MRI_num,'-0007*'], fullfile(save_directory,'Run_4'))
movefile(['rmeanf',MRI_num,'-0008*'], fullfile(save_directory,'Run_5'))

movefile(['wrmeanf',MRI_num,'-0004*'], fullfile(save_directory,'Run_1'))
movefile(['wrmeanf',MRI_num,'-0005*'], fullfile(save_directory,'Run_2'))
movefile(['wrmeanf',MRI_num,'-0006*'], fullfile(save_directory,'Run_3'))
movefile(['wrmeanf',MRI_num,'-0007*'], fullfile(save_directory,'Run_4'))
movefile(['wrmeanf',MRI_num,'-0008*'], fullfile(save_directory,'Run_5'))

movefile(['swrmeanf',MRI_num,'-0004*'], fullfile(save_directory,'Run_1'))
movefile(['swrmeanf',MRI_num,'-0005*'], fullfile(save_directory,'Run_2'))
movefile(['swrmeanf',MRI_num,'-0006*'], fullfile(save_directory,'Run_3'))
movefile(['swrmeanf',MRI_num,'-0007*'], fullfile(save_directory,'Run_4'))
movefile(['swrmeanf',MRI_num,'-0008*'], fullfile(save_directory,'Run_5'))
%}

%Move smoothed files
movefile(['swrrf',MRI_num,'-0005*'], fullfile(save_directory,'Run_1'))
movefile(['swrrf',MRI_num,'-0006*'], fullfile(save_directory,'Run_2'))
movefile(['swrrf',MRI_num,'-0007*'], fullfile(save_directory,'Run_3'))
movefile(['swrrf',MRI_num,'-0008*'], fullfile(save_directory,'Run_4'))
movefile(['swrrf',MRI_num,'-0009*'], fullfile(save_directory,'Run_5'))

%Move rp text files
movefile(['rp_f',MRI_num,'-0005*'], fullfile(save_directory,'Run_1'))
movefile(['rp_f',MRI_num,'-0006*'], fullfile(save_directory,'Run_2'))
movefile(['rp_f',MRI_num,'-0007*'], fullfile(save_directory,'Run_3'))
movefile(['rp_f',MRI_num,'-0008*'], fullfile(save_directory,'Run_4'))
movefile(['rp_f',MRI_num,'-0009*'], fullfile(save_directory,'Run_5'))

%Move mean nii files
movefile(['meanf',MRI_num,'-0005*'], fullfile(save_directory,'Run_1'))
movefile(['meanf',MRI_num,'-0006*'], fullfile(save_directory,'Run_2'))
movefile(['meanf',MRI_num,'-0007*'], fullfile(save_directory,'Run_3'))
movefile(['meanf',MRI_num,'-0008*'], fullfile(save_directory,'Run_4'))
movefile(['meanf',MRI_num,'-0009*'], fullfile(save_directory,'Run_5'))

movefile(['rmeanf',MRI_num,'-0005*'], fullfile(save_directory,'Run_1'))
movefile(['rmeanf',MRI_num,'-0006*'], fullfile(save_directory,'Run_2'))
movefile(['rmeanf',MRI_num,'-0007*'], fullfile(save_directory,'Run_3'))
movefile(['rmeanf',MRI_num,'-0008*'], fullfile(save_directory,'Run_4'))
movefile(['rmeanf',MRI_num,'-0009*'], fullfile(save_directory,'Run_5'))

movefile(['wrmeanf',MRI_num,'-0005*'], fullfile(save_directory,'Run_1'))
movefile(['wrmeanf',MRI_num,'-0006*'], fullfile(save_directory,'Run_2'))
movefile(['wrmeanf',MRI_num,'-0007*'], fullfile(save_directory,'Run_3'))
movefile(['wrmeanf',MRI_num,'-0008*'], fullfile(save_directory,'Run_4'))
movefile(['wrmeanf',MRI_num,'-0009*'], fullfile(save_directory,'Run_5'))

movefile(['swrmeanf',MRI_num,'-0005*'], fullfile(save_directory,'Run_1'))
movefile(['swrmeanf',MRI_num,'-0006*'], fullfile(save_directory,'Run_2'))
movefile(['swrmeanf',MRI_num,'-0007*'], fullfile(save_directory,'Run_3'))
movefile(['swrmeanf',MRI_num,'-0008*'], fullfile(save_directory,'Run_4'))
movefile(['swrmeanf',MRI_num,'-0009*'], fullfile(save_directory,'Run_5'))

%{
%Move smoothed files
movefile(['swrrf',MRI_num,'-0004*'], fullfile(save_directory,'Run_1'))
movefile(['swrrf',MRI_num,'-0005*'], fullfile(save_directory,'Run_2'))
movefile(['swrrf',MRI_num,'-0007*'], fullfile(save_directory,'Run_3'))
movefile(['swrrf',MRI_num,'-0008*'], fullfile(save_directory,'Run_4'))
movefile(['swrrf',MRI_num,'-0009*'], fullfile(save_directory,'Run_5'))

%Move rp text files
movefile(['rp_f',MRI_num,'-0004*'], fullfile(save_directory,'Run_1'))
movefile(['rp_f',MRI_num,'-0005*'], fullfile(save_directory,'Run_2'))
movefile(['rp_f',MRI_num,'-0007*'], fullfile(save_directory,'Run_3'))
movefile(['rp_f',MRI_num,'-0008*'], fullfile(save_directory,'Run_4'))
movefile(['rp_f',MRI_num,'-0009*'], fullfile(save_directory,'Run_5'))

%Move mean nii files
movefile(['meanf',MRI_num,'-0004*'], fullfile(save_directory,'Run_1'))
movefile(['meanf',MRI_num,'-0005*'], fullfile(save_directory,'Run_2'))
movefile(['meanf',MRI_num,'-0007*'], fullfile(save_directory,'Run_3'))
movefile(['meanf',MRI_num,'-0008*'], fullfile(save_directory,'Run_4'))
movefile(['meanf',MRI_num,'-0009*'], fullfile(save_directory,'Run_5'))

movefile(['rmeanf',MRI_num,'-0004*'], fullfile(save_directory,'Run_1'))
movefile(['rmeanf',MRI_num,'-0005*'], fullfile(save_directory,'Run_2'))
movefile(['rmeanf',MRI_num,'-0007*'], fullfile(save_directory,'Run_3'))
movefile(['rmeanf',MRI_num,'-0008*'], fullfile(save_directory,'Run_4'))
movefile(['rmeanf',MRI_num,'-0009*'], fullfile(save_directory,'Run_5'))

movefile(['wrmeanf',MRI_num,'-0004*'], fullfile(save_directory,'Run_1'))
movefile(['wrmeanf',MRI_num,'-0005*'], fullfile(save_directory,'Run_2'))
movefile(['wrmeanf',MRI_num,'-0007*'], fullfile(save_directory,'Run_3'))
movefile(['wrmeanf',MRI_num,'-0008*'], fullfile(save_directory,'Run_4'))
movefile(['wrmeanf',MRI_num,'-0009*'], fullfile(save_directory,'Run_5'))

movefile(['swrmeanf',MRI_num,'-0004*'], fullfile(save_directory,'Run_1'))
movefile(['swrmeanf',MRI_num,'-0005*'], fullfile(save_directory,'Run_2'))
movefile(['swrmeanf',MRI_num,'-0007*'], fullfile(save_directory,'Run_3'))
movefile(['swrmeanf',MRI_num,'-0008*'], fullfile(save_directory,'Run_4'))
movefile(['swrmeanf',MRI_num,'-0009*'], fullfile(save_directory,'Run_5'))
%}

%{
%Move smoothed files
movefile(['swrrf',MRI_num,'-0006*'], fullfile(save_directory,'Run_1'))
movefile(['swrrf',MRI_num,'-0007*'], fullfile(save_directory,'Run_2'))
movefile(['swrrf',MRI_num,'-0008*'], fullfile(save_directory,'Run_3'))
movefile(['swrrf',MRI_num,'-0009*'], fullfile(save_directory,'Run_4'))
movefile(['swrrf',MRI_num,'-0010*'], fullfile(save_directory,'Run_5'))

%Move rp text files
movefile(['rp_f',MRI_num,'-0006*'], fullfile(save_directory,'Run_1'))
movefile(['rp_f',MRI_num,'-0007*'], fullfile(save_directory,'Run_2'))
movefile(['rp_f',MRI_num,'-0008*'], fullfile(save_directory,'Run_3'))
movefile(['rp_f',MRI_num,'-0009*'], fullfile(save_directory,'Run_4'))
movefile(['rp_f',MRI_num,'-0010*'], fullfile(save_directory,'Run_5'))

%Move mean nii files
movefile(['meanf',MRI_num,'-0006*'], fullfile(save_directory,'Run_1'))
movefile(['meanf',MRI_num,'-0007*'], fullfile(save_directory,'Run_2'))
movefile(['meanf',MRI_num,'-0008*'], fullfile(save_directory,'Run_3'))
movefile(['meanf',MRI_num,'-0009*'], fullfile(save_directory,'Run_4'))
movefile(['meanf',MRI_num,'-0010*'], fullfile(save_directory,'Run_5'))

movefile(['rmeanf',MRI_num,'-0006*'], fullfile(save_directory,'Run_1'))
movefile(['rmeanf',MRI_num,'-0007*'], fullfile(save_directory,'Run_2'))
movefile(['rmeanf',MRI_num,'-0008*'], fullfile(save_directory,'Run_3'))
movefile(['rmeanf',MRI_num,'-0009*'], fullfile(save_directory,'Run_4'))
movefile(['rmeanf',MRI_num,'-0010*'], fullfile(save_directory,'Run_5'))

movefile(['wrmeanf',MRI_num,'-0006*'], fullfile(save_directory,'Run_1'))
movefile(['wrmeanf',MRI_num,'-0007*'], fullfile(save_directory,'Run_2'))
movefile(['wrmeanf',MRI_num,'-0008*'], fullfile(save_directory,'Run_3'))
movefile(['wrmeanf',MRI_num,'-0009*'], fullfile(save_directory,'Run_4'))
movefile(['wrmeanf',MRI_num,'-0010*'], fullfile(save_directory,'Run_5'))

movefile(['swrmeanf',MRI_num,'-0006*'], fullfile(save_directory,'Run_1'))
movefile(['swrmeanf',MRI_num,'-0007*'], fullfile(save_directory,'Run_2'))
movefile(['swrmeanf',MRI_num,'-0008*'], fullfile(save_directory,'Run_3'))
movefile(['swrmeanf',MRI_num,'-0009*'], fullfile(save_directory,'Run_4'))
movefile(['swrmeanf',MRI_num,'-0010*'], fullfile(save_directory,'Run_5'))
%}

toc

%% Clear intermediate file types

cd(data_directory)
delete rf* rrf* wrrf*
