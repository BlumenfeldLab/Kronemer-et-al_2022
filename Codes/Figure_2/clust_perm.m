%% EEG Channel Based Permutation Analysis - ROI Timecourse

%The purpose of this code is to take voxels from an ROI, average over those
%voxels (thus collapsing the space dimension), and find statically
%signficiant changes in time via the cluster based permutation method. This
%is completed both by CP vs CnP and CP-CnP vs baseline (prestimulus period)

%Written by: Sharif I. Kronemer
%Date: 3/10/2021
%Modified: Thomas Xin 6/28/2021

clear

%% Run Location

%Select run location
prompt_1 = 'Running code local or server [l, s]: ';
run_location = input(prompt_1,'s');

use_block = false ; 
len_block = 50; 
prompt_2 = "Running EEG or RNS analyis [eeg, rns]: "; 
data_type = input(prompt_2,'s'); 

%% Directories and Variable Names

%Variable name
folder_name = 'cent_quad_threshold_1_15s';

if isequal(run_location, 's')
    
    %Add behavioral analysis path
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/Behavioral Analysis')
    
    %Add paths
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Permutation Statistical Analysis/Permutation functions')
    addpath('/mnt/Data5/Sharif Data F/HNCT fMRI Study/fMRI Analysis/Analysis Codes/Time courses/Supplementary functions/fMRI Analysis/Code for Percent Change Maps/subRoutines')
    addpath('/mnt/Data8/HNCT_NRP_Study/HNCT No Report Paradigm/Analysis/Analysis Code/MRI Analysis/Permutation Statistical Analysis')

    %Data directory
    data_dir = '/mnt/Data7/RNS Study/Analysis/Analysis Code/Temporal Analysis';
    
    %Save directory
    save_dir = '/mnt/Data7/RNS Study/Analysis/Group Analysis/Voltage ERP';
         
elseif isequal(run_location, 'l')
       
end

%% Load group voltage data

cd(data_dir)
data_name = ['group_' data_type '.mat']; 
group_data = load(data_name); 

if strcmp(data_type,'rns') 
    
    place_holder = zeros(1,1001,353,7); 
    group_data.group_eeg_cp = myReshape(cat(1,nanmean(group_data.group_eeg_cp,1),place_holder)); 
    group_data.group_eeg_cnp = myReshape(cat(1,nanmean(group_data.group_eeg_cnp,1),place_holder)); 

else 
    
    group_data.group_eeg_cp = myReshape(group_data.group_eeg_cp); 
    group_data.group_eeg_cnp = myReshape(group_data.group_eeg_cnp); 
    
end 


if ~use_block
    
    block_name = []; 
    [group_data.group_eeg_cp,group_data.group_eeg_cnp] = selectTrials(group_data.group_eeg_cp,group_data.group_eeg_cnp);

else 
    
    block_name = 'block'; 
    max_len = min(size(group_data.group_eeg_cp,3),size(group_data.group_eeg_cnp,3)); 
    group_data.group_eeg_cp = getBlocks(group_data.group_eeg_cp,len_block,max_len); 
    group_data.group_eeg_cnp = getBlocks(group_data.group_eeg_cnp,len_block,max_len); 
    
end 

%% Parameters

%Number of permutations
num_permutations = 10000;

%Are samples dependent (default is true)
dependent_samples = 'true';

%Define alpha threshold
p_threshold = 0.05;

%Two-sided
two_sided = 'true';

%% Main and Subtraction Epoch Subtraction and Baseline 

disp('Baselining data from prestim period')

%Specify baseline period 
baseline_window = 251:500;

%CP minus CnP

%Subtract main data from subtract data
group_pc_data = group_data.group_eeg_cp - group_data.group_eeg_cnp;

%Calculate baseline values [voxel x time x subjects]
group_pc_baseline = squeeze(nanmean(group_pc_data(:,baseline_window,:),2));

%Convert baseline [channel x time x subjects]
group_CP_minus_CnP_baseline = permute(repmat(group_pc_baseline,[1,1,1001]),[1,3,2]);

%Subtract baseline from main data
group_CP_minus_CnP_baselined_data = group_pc_data - group_CP_minus_CnP_baseline;

%Mean voltage over subjects
CP_minus_CnP_mean_voltage = nanmean(group_CP_minus_CnP_baselined_data,3);

%CP

%Calculate baseline values [channel x subjects]
group_pc_baseline = squeeze(nanmean(group_data.group_eeg_cp(:,baseline_window,:),2));

%Convert baseline [channel x time x subjects]
group_pc_baseline = permute(repmat(group_pc_baseline,[1,1,1001]),[1,3,2]);

group_CP_baseline = group_pc_baseline; 

%Subtract baseline from main data
group_CP_baselined_data = group_data.group_eeg_cp - group_pc_baseline;

%Mean voltage over subjects
CP_mean_voltage = nanmean(group_CP_baselined_data,3);

%CnP
    
%Calculate baseline values [voxel x time x subjects]
group_pc_baseline = squeeze(nanmean(group_data.group_eeg_cnp(:,baseline_window,:),2));

%Convert baseline [channel x time x subjects]
group_pc_baseline = permute(repmat(group_pc_baseline,[1,1,1001]),[1,3,2]);

%Subtract baseline from main data
group_CnP_baselined_data = group_data.group_eeg_cnp - group_pc_baseline;

%Mean voltage over subjects
CnP_mean_voltage = nanmean(group_CnP_baselined_data,3);

%% Run Permutation Tests

if strcmp(data_type,'eeg') 
    channel_list = [6,10,13,17];
    channel_name = {'Fz','Cz','Pz','Oz'}; 
    
elseif strcmp(data_type,'rns') 
    channel_list = [1]; 
    channel_name = {'RNS'};
    
end 

%Loop over channels
for chan = 1:length(channel_list)
    
    tic
    
    %Define channel
    current_channel = channel_list(chan);
    current_channel_name = channel_name{chan};

    %CP vs CnP Testing
    disp(['Running Permutation Tests - ',num2str(current_channel)])
    
    %CP vs CnP Testing
    [clusters, pval, t_sums, permutation_distribution] = permutest_TimeCourses(squeeze(group_CP_baselined_data(current_channel,:,:)), squeeze(group_CnP_baselined_data(current_channel,:,:)), dependent_samples, ...
        p_threshold, num_permutations, two_sided);

    %Find significant clusters pvalue < 0.05
    sig_clust = find(pval < 0.05);

    %Find the significant time points
    sig_time_pts = sort([clusters{sig_clust}]);

    %Save output
    cd(save_dir)
    save(['timecourse_cluster_CP_vs_CnP_channel_',num2str(current_channel),'_',current_channel_name,'_',num2str(num_permutations),'perm.mat'],...
        'clusters','pval','t_sums','permutation_distribution','sig_clust','sig_time_pts');

    %CP vs Baseline Testing
    [clusters, pval, t_sums, permutation_distribution] = permutest_TimeCourses(...
        squeeze(group_CP_baselined_data(current_channel,:,:)),...
        squeeze(group_CP_baseline(current_channel,:,:)), dependent_samples, ...
        p_threshold, num_permutations, two_sided);

    %Find significant clusters pvalue < 0.05
    sig_clust = find(pval < 0.05);

    %Find the significant time points
    sig_time_pts = sort([clusters{sig_clust}]);

    %Save output
    cd(save_dir)
    save(['timecourse_cluster_CP_minus_CnP_channel_',num2str(current_channel),'_',current_channel_name,'_',num2str(num_permutations),'perm.mat'],...
        'clusters','pval','t_sums','permutation_distribution','sig_clust','sig_time_pts');
  
    toc
    
end

%% Plot Timecourse ERPs and Signifance Bars

%Set renderer to painter
set(0, 'DefaultFigureRenderer', 'painters');

%ERP Times
N100_on = [75,75];
N100_off =[125,125];

VAN_on = [175,175];
VAN_off = [225,225];

P2_on = [275,275];
P2_off = [325,325];

P3_on = [350,350];
P3_off = [650,650];

%Start and End Times
stim_time = 501;
start_time = -500;
end_time = 1500;

%Axis limits/voltage scale
max_voltage = 15;
min_voltage = -10;

%Loop over channels
for chan = 1:length(channel_list)

    %Define current channel
    current_channel = channel_list(chan);
    current_channel_name = channel_name{chan};

    figure
    hold on

    %Figure labels
    title(['ERPs  - Channel ', num2str(current_channel), ' ',current_channel_name], 'Interpreter', 'none')
    xlabel('Time (ms)')
    ylabel('Voltage (uV)')
    ylim([min_voltage max_voltage])
    xlim([start_time end_time])

    %Define timevector
    timevector_plot = [start_time:4:end_time];
    timevector_data = [501-abs(start_time)/4:501+ceil(end_time/4)];

    %Reference lines
    plot([0,0],[min_voltage,max_voltage],'k')
    plot([start_time, end_time],[0,0],'k')
    
    %N100
    plot(N100_on,[min_voltage,max_voltage],':k')
    plot(N100_off,[min_voltage,max_voltage],':k')
    
    %VAN
    plot(VAN_on,[min_voltage,max_voltage],':k')
    plot(VAN_off,[min_voltage,max_voltage],':k')
    
    %P2
    plot(P2_on,[min_voltage,max_voltage],':k')
    plot(P2_off,[min_voltage,max_voltage],':k')
    
    %P3
    plot(P3_on,[min_voltage,max_voltage],':k')
    plot(P3_off,[min_voltage,max_voltage],':k')
    
    %%True class
    CP = plot(timevector_plot, CP_mean_voltage(current_channel,timevector_data),'b')
    CnP = plot(timevector_plot, CnP_mean_voltage(current_channel,timevector_data),'r')
    CP_minus_CnP = plot(timevector_plot, CP_minus_CnP_mean_voltage(current_channel,timevector_data),'g')
    
    %Plot sig bars 
    bar_types = {'CP_vs_CnP','CP_minus_CnP'};%'CP_vs_PP','CnP_vs_PnP'};
    bar_colors = {'m','k'};
    bar_yaxis = [12,10];
    
    %Loop over bar types
    for bar = 1:length(bar_types)

        %Define current bar
        current_bar = bar_types{bar};

        %Define current bar color
        current_color = bar_colors{bar};

        %Define current bar axis
        current_yaxis = bar_yaxis(bar);

        %Load bars
        load(fullfile(save_dir,['timecourse_cluster_',current_bar,'_channel_',num2str(current_channel),'_',current_channel_name,'_',num2str(num_permutations),'perm.mat']))

        %If there are significant times points
        if not(isempty(sig_time_pts))

            %Find significant break points
            sig_breaks = find(diff(sig_time_pts) > 1);
            
            %If there are break points
            if not(isempty(sig_breaks))
            
                %Loop over discontinuous points
                for breaks = 1:length(sig_breaks)

                    %Current break point
                    break_point = sig_breaks(breaks);

                    %First break
                    if breaks == 1 

                        plot((sig_time_pts(1:break_point)-stim_time)*4,ones(1,length(sig_time_pts(1:break_point)))*current_yaxis, 'Color', current_color, 'LineWidth',2)

                        %If only one break point also plot to the end of
                        %sig times
                        if length(sig_breaks) == 1

                            %Between last break and end
                            plot(([sig_time_pts(break_point+1:length(sig_time_pts))]-stim_time)*4,ones(1,length(sig_time_pts(break_point+1:length(sig_time_pts))))*current_yaxis, 'Color', current_color, 'LineWidth',2)

                        end

                    %Last break
                    elseif breaks == length(sig_breaks)

                        %Between penultimate and final break
                        plot(([sig_time_pts(sig_breaks(breaks-1)+1):sig_time_pts(break_point)]-stim_time)*4,...
                            ones(1,length((sig_time_pts(sig_breaks(breaks-1)+1):sig_time_pts(break_point))))*current_yaxis, 'Color', current_color, 'LineWidth',2)

                        %Between last break and end
                        plot(([sig_time_pts(break_point+1:length(sig_time_pts))]-stim_time)*4,ones(1,length(sig_time_pts(break_point+1:length(sig_time_pts))))*current_yaxis, 'Color', current_color, 'LineWidth',2)

                    %Middle breaks
                    else

                        plot(([sig_time_pts(sig_breaks(breaks-1)+1):sig_time_pts(break_point)]-stim_time)*4,...
                            ones(1,length((sig_time_pts(sig_breaks(breaks-1)+1):sig_time_pts(break_point))))*current_yaxis, 'Color', current_color, 'LineWidth',2)

                    end

                end
            
            %If there are no breaks in the sig time points    
            else
                
                %Plot sig times
                plot(([sig_time_pts]-stim_time)*4,ones(1,length(sig_time_pts))*current_yaxis, 'Color', current_color, 'LineWidth',2)

            end

        end

    end

    %Define legend
    legend([CP, CnP, CP_minus_CnP], ['CP']', ['CnP'], ['CP minus CnP'], 'Location', 'best')

    cd(save_dir)                      
    savefig(['ERPs_CP_CnP_channel_',num2str(current_channel),'_',current_channel_name,'_',block_name,'.fig'])

    close

end

%% Define functions

function output = getBlocks(input1,lenBlock,max_len)
    len1 = size(input1,3); 
    
    niter = ceil(max_len/lenBlock); 
    
    output = [] ; 
    for ii = 1:niter
        if ii ~= max(niter)
            one_block = input1(:,:,(ii-1)*lenBlock+1:ii*lenBlock); 
        else 
            one_block = input1(:,:,(ii-1)*lenBlock+1:end); 
        end 
        output = cat(3,output,nanmean(one_block,3)); 
    end 
    
end 

function [output1,output2] = selectTrials(input1,input2)
    len1 = size(input1,3); 
    len2 = size(input2,3); 
    
    min_len = min(len1,len2); 
    max_len = max(len1,len2); 
        
    perm_idx1 = randperm(max_len,min_len) ; 
    perm_idx2 = randperm(min_len); 
    
    output1 = input1(:,:,perm_idx1); 
    output2 = input2(:,:,perm_idx2); 
end 

function output = myReshape(input)
[x,y,z,w] =size(input) ; 
    output = []; 
    for iz = 1:z
        for iw = 1:w
            one_arr  = squeeze(input(:,:,iz,iw)); 
            if any(all(isnan(one_arr),2))
                continue 
            end 
            output = cat(3,output,one_arr);
        end 
    end 
end 