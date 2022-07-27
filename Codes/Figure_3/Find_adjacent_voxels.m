%% Finding "grid" voxels and "adjacent" voxels

%The purpose of this code is to test approaches of identifying voxels of
%interest among a subset of brain voxels. 

%Critically this code creates the adjacency voxel matrix that will be used
%in permutation testing. This matrix is built from the 91x109x91

%Written by: Sharif I. Kronemer
%Date: 6/2/2020

clear

%% Directories

%Whole brain grey matter masked voxels
ROI_dir = '/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Subject Analysis/191AM/MRI Analysis/Extracted voxel data/Run 1 Movie/axialSlices'; %MNI brain
%ROI_dir = 'Z:\Sharif Data F\HNCT fMRI Study\fMRI Analysis\Subject Analysis\191AM/MRI Analysis\Extracted voxel data\Run 1 Movie\axialSlices'; %MNI brain

%Save directory
save_dir = '/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Permutation Statistical Analysis'; 
%save_dir = 'S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis/Analysis Code\MRI Analysis\Permutation Statistical Analysis';

%% Load Grey Masked Brain

%Load voxelInfo for that ROI
cd(ROI_dir)
load(fullfile(ROI_dir,'voxelInfo.mat'));

%Define voxel coordinates     
vXYZ = xyz.vXYZ; 

%Image dimensions - Default 91x109x91
sizeMat = [91,109,91];

%Define index of voxel XYZ positions for image of size 91x109x91
alloutInd = mat2ind(vXYZ,sizeMat);

%Adjacent voxels - all possible adjacent voxel positions
base = [+1 +1 0; +1 -1 0; +1 +1 +1; +1 0 +1; +1 -1 +1; +1 +1 -1; ...
+1 0 -1; +1 -1 -1; +1 0 0; 0 +1 0; 0 -1 0; 0 +1 +1; 0 0 +1; 0 -1 +1; ...
0 +1 -1; 0 0 -1; 0 -1 -1; -1 +1 0; -1 -1 0; -1 +1 +1; -1 0 +1; -1 -1 +1; ...
-1 +1 -1; -1 0 -1; -1 -1 -1; -1 0 0];

%Define max and min dimensions of grey matter masked img
x_min = 11;
x_max = 82;

y_min = 12;
y_max = 100;

z_min = 7;
z_max = 78;

%% Creating Full Box Matrix 10-82 (x) 11-100 (y) 6-78 (z)
%Note: Makes a 3D rectangle of positions within the bounds of the grey
%matter masked brain.

%{
%X-axis loop counter
count = 0;

%Loop over x-axis 
for x = 10:82
    
    disp(['Running ', num2str(x)])
    
    %Loop over y-axis
    for y = 11:100
        
        %Loop over z-axiz
        for z = 6:78
            
            %Add to counter
            count = count +1;
            
            %Store X,Y,Z combination
            box_XYZ(:,count) = [x;y;z];
                        
        end
    end
end

%Sort dimensions according to z-axis
[temp, order] = sort(box_XYZ(3,:));
box_XYZ = box_XYZ(:, order);
%}

%% Select Grid Voxels

disp('Determine Grid Voxels')

%X-axis loop counter
count = 0;

%Loop over x-axis 
for x = x_min:3:x_max
    
    %disp(['Running ', num2str(x)])
    
    %Loop over y-axis
    for y = y_min:3:y_max
        
        %Loop over z-axiz
        for z = z_min:3:z_max
    
            %Add to counter
            count = count +1;
            
            %Store X,Y,Z combination
            grid_XYZ(:,count) = [x;y;z];
            
        end
        
    end
    
end

%Sort dimensions according to z-axis
[temp, order] = sort (grid_XYZ(3,:)); %Get order on z-axis
grid_XYZ = grid_XYZ(:, order); %Sort according to order

%% Find Grid Voxels Among Brain Voxels in Brain Mask

disp('Find Grid Pts Among Brain Voxels')

tic

%Convert grid voxel array to cells
cell_grid_XYZ = num2cell(grid_XYZ,1);

%Convert grey matter mask voxel array to cells
brain_XYZ = num2cell(vXYZ,1);

%Initialize variable
brain_grid_voxels = [];

%Loop over the number of grid voxels
for grid_idx = 1:size(cell_grid_XYZ,2)
    
    disp(['Running ', num2str(grid_idx)])
   
    %Loop over all brain voxels
    for brain_idx = 1:size(brain_XYZ,2)
        
        %Check if brain voxel XYZ matches grid voxel XYZ
        if isequal(cell_grid_XYZ(grid_idx), brain_XYZ(brain_idx))
            
            %Store brain voxel index
            brain_grid_voxels(end+1) = brain_idx;
            
            %Skip to next grid voxel once corresponding brain voxel is found
            break
        
        end
        
    end
    
end

%Sort the grid voxels for organization
brain_grid_voxels = sort(brain_grid_voxels);

toc

%% Find adjacent voxels to grid voxel positions found among brain voxels

disp('Find brain-grid voxel adjacent voxels')

%Setup variable
all_adjacent_voxels_array = [];
all_adjacent_voxels_cell = {};

tic

%Loop over brain-grid voxels
for bg_voxel = 1:length(brain_grid_voxels)
    
        disp(['Running voxel ', num2str(bg_voxel)])
        
        %Select current voxel
        test_xyz = vXYZ(:,brain_grid_voxels(bg_voxel));
        
        %Adjacent voxel count to bg_voxel
        adjacent_voxel_count = 0;
        
        %Reset voxel list
        voxel_list = [];
        
        %Loop over possible adjacent locations
        for loc = 1:length(base)
            
            %Define an adjacent voxel location
            adjacent_voxel = sum([base(loc,:)', test_xyz], 2);
            
            %Check if the adjacent voxel exists in brain mask voxels by comparing the XYZ positions
            if ~isempty(find(adjacent_voxel(1) == vXYZ(1,:) & adjacent_voxel(2) == vXYZ(2,:) & adjacent_voxel(3) == vXYZ(3,:)))
                
                %Add to adjacent voxel count
                adjacent_voxel_count = adjacent_voxel_count +1;
                
                %Find adjacent voxel index in 91x109x91 matrix and store
                voxel_list(adjacent_voxel_count) = mat2ind(adjacent_voxel,sizeMat); 
                
            end
            
        end
        
        %Add cell of adjacent voxels in cell
        all_adjacent_voxels_cell{bg_voxel} = voxel_list';
        
        %Add array of adjacent voxels in matrix
        all_adjacent_voxels_array = [all_adjacent_voxels_array; voxel_list'];

end

toc

%% Save grid and adjacent voxels matrices

%Save
cd(save_dir)
save('brain_grid_voxel_info.mat', 'brain_grid_voxels', 'all_adjacent_voxels_array', 'all_adjacent_voxels_cell', '-v7.3')

%% Plot Point Cloud

%Setup figure
figure
hold on

title('Voxel Point Cloud')
xlabel('X')
ylabel('Y')
zlabel('Z')

%Loop over voxels to plot
for voxel = 20000:25000
    
   %Check if grid voxels
   if ismember(voxel, brain_grid_voxels)
       
        plot3(vXYZ(1,voxel),vXYZ(2,voxel), vXYZ(3,voxel), 'o', 'Color', 'r', 'MarkerFaceColor', 'r')
        
   %Check if adjacent voxels
   elseif ismember(voxel, all_adjacent_voxels_array)
       
       plot3(vXYZ(1,voxel),vXYZ(2,voxel), vXYZ(3,voxel), 'o', 'Color', 'k', 'MarkerFaceColor', 'y')
   
   %All other voxels
   else
       
       plot3(vXYZ(1,voxel),vXYZ(2,voxel), vXYZ(3,voxel), 'o', 'Color', 'k')
       
   end 
   
end
