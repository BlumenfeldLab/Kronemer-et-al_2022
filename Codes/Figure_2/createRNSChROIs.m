%% Use MarsBar to make RNS Channel ROIs

%This code will compute the channel location from the contact coordinates
%and generate an ROI at this position. This code will also combine ROIs
%based on above or below the 3uv threshold post-stimulus.

%Written by: Sharif I. Kronemer Edited by Thomas Xin 
%Date: 5\22\2021

clear 

%% Add paths

addpath(['V:\Sharif Data F\HNCT fMRI Study\fMRI Analysis\Analysis Codes\marsbar-0.44']);
addpath('V:\Sharif Data F\HNCT fMRI Study\spm12_radiological')

%% big brain template file 

imgname = fullfile('R:\RNS Study\RNS Locations\Combined_ROI','bigbrain-temp','PD25-MPRAGET1w-template-300um.nii'); 

%% create ROIs for each subject

%Subject ID
all_subjects = {'471MH','489RD','491GS','490KB','535BP','536JD','451NH'};
all_flatten_dims = {'3D'};

channel_table = []; 
comb_save_path = 'R:\RNS Study\RNS Locations\Combined_ROI\RNS Channels\3D-auto\'; 

%Loop over subjects
for isub = 1:length(all_subjects)
    
    for idim = 1:length(all_flatten_dims)
        
        one_table = createROIbySubject(all_subjects{isub},all_flatten_dims{idim},comb_save_path,mars_space(imgname)); 
        channel_table = [channel_table; one_table];
     
    end
    
end 

%% sort output tables 
channel_table = sortrows(channel_table,'Subject');

%% Define Funtion
function channel_table = createROIbySubject(ID,flatten_dim,comb_save_path,image_space)

    %% create table for manual change of coordiantes in case of overlap
    change_table = array2table(zeros(5,4));
    change_table.Properties.VariableNames = {'Subject','Channel','Hemi','Change'}; 
    change_table.Subject = {'490KB','471MH','490KB','490KB','535BP'}'; 
    change_table.Hemi = {'L','L','L','R','L'}'; 
    change_table.Channel = [2,1, 1, 1,2]'; 
    change_table.Change = [0,-3,-3,3,3]';
    change_table.Dim = {'3D','coronal','coronal','coronal','coronal'}';

    group_names = {'above','below'}; 
    
    %read table for above and below threshold differences 
    
    channel_table = array2table(zeros(4,5));

    channel_table.Properties.VariableNames = {'Subject','Channel','X','Y','Z'}; 
    channel_table.Subject = repmat({ID},4,1); 
    channel_table.Channel = {'1L','2L','1R','2R'}'; 
    nh_ch_coord = {[12, -16, 8],[15,-19,12]};
    
    MNI_savename = {'sagittal','coronal','axial'};
    
    MNI_coord.sagittal = 8;
    MNI_coord.coronal = -20; 
    MNI_coord.axial = 1;
   
    % read rns table in workspace 
    read_RNS_contact_coordinates;

    addpath('\\172.23.254.106\Data7\RNS Study\RNS Locations')
  
    RNS_location_filepath = 'R:\RNS Study\RNS Locations\Blumenfeld Lab RNS_coordinates.xlsx';
    
    %% Set general optionstion_
    %Flatten coorindate
    if strcmp(flatten_dim,'sagittal') 
            dim = 1; 
    elseif strcmp(flatten_dim,'coronal') 
            dim = 2; 
    elseif  strcmp(flatten_dim,'axial')
            dim = 3; 
    else
            dim = NaN; 
    end 
    
    %Out put directory

     outDir = fullfile('R:\RNS Study\RNS Locations',ID,'Channel ROI',flatten_dim); %Data7

    mkdir(outDir)
    cd(outDir)

    %Size (mm)
    sphereRadius = 1;
    boxwidths = [1 1 1];%[4 4 4]

    %% Create channel ROI on right and left hemisphere

    %Hemisphere
    if ~isequal(ID,'451NH')
        hemispheres = {'R','L'};
    else 
        hemispheres = {'N'}; 
    end 

    %Loop over hemisphere
    for side = 1:length(hemispheres)

        %Select hemisphere
        hemi = hemispheres{side};
        
        if isequal(ID,'451NH')

            contacts = readtable(RNS_location_filepath,'Sheet','RNS Channel Location',...
                'Range','B28:G29','ReadVariableNames',true);
                    %Find channel 1 location - midpoint of 3D contact locations
                    
            for ich = 1:2
                
                    channel_coord{ich} = contacts{ich,3:5}; 
            end 
        else
            %Select the contact locations
            if isequal(hemi,'R')
                contacts = RNS_location_table(strcmpi(RNS_location_table.subject,ID),2:4); 


            elseif isequal(hemi,'L')
                contacts = RNS_location_table(strcmpi(RNS_location_table.subject,ID),5:7); 


            end
        
            ch_done = [] ; 
            channel_coord = {}; 
            for ich = 1:size(contacts,1)
                if ~any(ch_done==ich) 
                channel_coord = [channel_coord  {(contacts{ich,:} + contacts{ich+1,:})/2}];
                ch_done = [ch_done ich+1];
                end 
            end 
           
        end 
        
        %Flatten on the axial dimension
        if ~isnan(dim)

            %Replace Z coordinate with constant coordinate
            channel_coord{1}(dim) = MNI_coord.(MNI_savename{dim});
            channel_coord{2}(dim) = MNI_coord.(MNI_savename{dim});

            %Save name
            save_name = ['2D' MNI_savename{dim}];

        else

            %Save name
            save_name = '3D';

        end
        
        %% Make ROIs and Nifiti Files

        fprintf('\n');
       
        %Loop over coordinates
        for channel = 1:length(channel_coord) 

            %Select current voxel
            thisCoord = channel_coord{channel}; 
            
            % change coordiantes for those with ovelap 
            if ismember(ID,change_table.Subject)
                change_table = change_table(strcmp(change_table.Dim,flatten_dim),:); 
                
                change = change_table(strcmp(change_table.Subject,ID) & change_table.Channel==channel & strcmp(change_table.Hemi,hemi),'Change'); 
                if isempty(change) 
                    thisCoord(1)= thisCoord(1);  
                else 
                thisCoord(1) = thisCoord(1) + change{1,1}; 
                end 
            end 

            %Define voxel coordinates
            roiLabel = sprintf('%i_%i_%i', thisCoord(1), thisCoord(2), thisCoord(3));
            
            %Generate spherical ROI
            ROI{side,channel} = maroi_sphere(struct('centre', round(thisCoord), 'radius', sphereRadius));
            %ROI{side,channel} = maroi_box(struct('centre', round(thisCoord), 'widths', boxwidths)); 

            %Define directory
            outName = fullfile(outDir, ['Channel ',num2str(channel),hemi,'_',save_name]);
             
            if ~isequal(ID,'451NH')
                try
                channel_table{strcmp(channel_table.Channel,[num2str(channel) hemi]),3:5}= thisCoord; 
                catch
                    disp('la')
                end
   
            end 
                        
            %Save MarsBaR ROI (.mat) file
            saveroi(ROI{side,channel}, [outName '_roi.mat']);

            %Save the Nifti (.nii) file using big brain templates
            mars_rois2img(ROI{side,channel}, [outName '.nii'], image_space,'i')

        end

    end
    
    %% combine and save ROIs  
    
    % iterate through above and below threshold groups 
    for ig = 1:2
        group_name = group_names{ig}; 
        one_table = group_table(strcmpi(group_table.Subject,ID)&strcmp(group_table.Group,group_name),:); 
        ch_idx = one_table.Channel; 
        if isempty(one_table)
            continue
        end 
        
        hemi_idx = cellfun(@(x) find(strcmp(x,hemispheres)), one_table.Hemisphere); 
        
        % combine ROIs
        func = 'r1'; 
        for i = 1:length(ch_idx)
            eval(sprintf('r%d = ROI{hemi_idx(%d),ch_idx(%d)};', i, i,i));
            if i ~=1 
                func = [func '*' sprintf('r%d',i)] ;
            end 
        end
        
        % get Marsbar setup 
        [Finter,Fgraph,CmdLine] = spm('FnUIsetup','Combine ROIs');
         eval(['o=' func ';']);

        %Save MarsBaR ROI (.mat) file
        combined_name = [ID '_' flatten_dim '_' group_name '_combined']; 
        saveroi(o, [comb_save_path combined_name '_roi.mat']);

        %Save the Nifti (.nii) file
         mars_rois2img(o, [comb_save_path combined_name '.nii'], image_space,'i')

    end
    
end 