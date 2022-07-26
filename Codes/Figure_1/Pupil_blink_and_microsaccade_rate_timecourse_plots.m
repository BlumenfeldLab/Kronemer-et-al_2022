%% Plot Pupil, Blink, and Microsaccade Timecourses

%This code will load the group Eyelink data and the group eyelink data
%permutation statistical test results and plot the eyelink timecourses 

%Written by: Sharif I. Kronemer
%Date: 5/31/2021
%Modified: 7/4/2021

clear

%% Prompt

%Predicted or Confirmed Trials
prompt_1 = 'Predicted (PP/PnP) or confirmed (CP/CnP) trials [pt, ct]: ';
trial_type = input(prompt_1, 's');

%Location relevant
prompt_2 = 'No-report relevant location [c, q, combine]: ';
relevant_location = input(prompt_2, 's');

%Opacity
prompt_3 = 'Opacity [threshold, opaque, blank]: ';
opacity = input(prompt_3,'s');

%Imaging modality
prompt_4 = 'Modality type [EEG, MRI, both]: ';
modality = input(prompt_4,'s');

%EEG
if isequal(modality,'EEG')
   
    modality_name = 'EEG';
    
%MRI
elseif isequal(modality,'MRI')
    
    modality_name = 'MRI';
    
%EEG and MRI
elseif isequal(modality,'both')
    
    modality_name = 'EEG_MRI';
    
end

%Select confidence score if predicted trial type
if isequal(trial_type,'pt')
    
    %Confidence threshold value (Note: class 1 uses the positive value and
    %class 2 uses the negative value)
    prompt_5 = 'Confidence threshold value [0,0.25,0.5,0.75,1]: ';
    confidence_score = input(prompt_5,'s'); 
    
    %Select relevant or irrelevant condition
    prompt_6 = 'Relevant or irrelvant condition [relevant, irrelevant, rel_irrel]: ';
    relevant_condition = input(prompt_6,'s');
       
end

%With trials with blink rejected
prompt_7 = 'Blink trials rejected [y,n]: ';
blink_trials_rejected = input(prompt_7,'s');

%Define relevant condition variable for report trials
if isequal(trial_type,'ct')
    
    relevant_condition = 'relevant';
 
end

%Center Relevant/Irrelevant
if isequal(relevant_location, 'c')
  
    if isequal(relevant_condition,'relevant')
        
        relevant_location = 'Center Relevant';
    
    elseif isequal(relevant_condition,'irrelevant')
        
        relevant_location = 'Center Irrelevant';
        
    elseif isequal(relevant_condition,'rel_irrel')
        
        relevant_location = 'Center Relevant and Irrelevant';
        
    end 
        
    rel_save_name = 'cent';

%Quadrant Relevant
elseif isequal(relevant_location, 'q')
    
    if isequal(relevant_condition,'relevant')
        
        relevant_location = 'Quadrant Relevant';
    
    elseif isequal(relevant_condition,'irrelevant')
        
        relevant_location = 'Quadrant Irrelevant';
    
    elseif isequal(relevant_condition,'rel_irrel')

        relevant_location = 'Quadrant Relevant and Irrelevant';
        
    end 
    
    rel_save_name = 'quad';
    
%Combine Center and Quadrant Relevant   
elseif isequal(relevant_location, 'combine')
        
    if isequal(relevant_condition,'relevant')
        
        relevant_location = 'Center and Quadrant Relevant';
    
    elseif isequal(relevant_condition,'irrelevant')
        
        relevant_location = 'Center and Quadrant Irrelevant';
    
    elseif isequal(relevant_condition,'rel_irrel')
        
        relevant_location = 'Center and Quadrant Relevant and Irrelevant';

    end 
    
    rel_save_name = 'cent_quad';
    
end

%% Directories 

%Eyelink directory
eyelink_dir = 'S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Group Data\EyeLink';

%Eyelink filename
if isequal(trial_type,'ct')
    
    if isequal(blink_trials_rejected,'y')

        eyelink_filename = ['Group_eyelink_CP_CnP_',rel_save_name,'_',opacity,'_',modality_name,'_data.mat'];
    
    else
        
        eyelink_filename = ['Group_eyelink_CP_CnP_',rel_save_name,'_',opacity,'_',modality_name,'_no_blink_rejections_data.mat'];
        
    end
    
elseif isequal(trial_type,'pt')
    
    if isequal(blink_trials_rejected,'y')
   
        eyelink_filename = ['Group_eyelink_',relevant_condition,'_PP_PnP_',rel_save_name,'_',opacity,'_',modality_name,'_score_thres_',num2str(confidence_score),'_data.mat'];

    else
        
        eyelink_filename = ['Group_eyelink_',relevant_condition,'_PP_PnP_',rel_save_name,'_',opacity,'_',modality_name,'_score_thres_',num2str(confidence_score),'_no_blink_rejections_data.mat'];
        
    end
    
end

%Statistical result directory
if isequal(trial_type,'ct')

    if isequal(blink_trials_rejected,'y')

        stat_dir = fullfile('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Group Analysis\EyeLink Analysis\Permutation Analysis',['CP_CnP_',rel_save_name,'_',opacity,'_',modality_name]);

    else

        stat_dir = fullfile('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Group Analysis\EyeLink Analysis\Permutation Analysis',...
            ['CP_CnP_',rel_save_name,'_',opacity,'_',modality_name,'_no_blink_rejections']);
       
    end
        
elseif isequal(trial_type,'pt')
    
    if isequal(blink_trials_rejected,'y')

        stat_dir = fullfile('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Group Analysis\EyeLink Analysis\Permutation Analysis',...
            ['PP_PnP_',relevant_condition,'_',rel_save_name,'_',opacity,'_',modality_name,'_confid_score_',num2str(confidence_score)]);

    else

        stat_dir = fullfile('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Group Analysis\EyeLink Analysis\Permutation Analysis',...
            ['PP_PnP_',relevant_condition,'_',rel_save_name,'_',opacity,'_',modality_name,'_confid_score_',num2str(confidence_score),'_no_blink_rejections']);

    end
        
end

%Save directory

%Filename
if isequal(trial_type,'ct')
    
    if isequal(blink_trials_rejected,'y')

        save_filename = ['CP_CnP_',rel_save_name,'_',opacity,'_',modality_name];
    
    else
        
        save_filename = ['CP_CnP_',rel_save_name,'_',opacity,'_',modality_name,'_no_blink_rejections'];

    end
    
elseif isequal(trial_type,'pt')
    
    if isequal(blink_trials_rejected,'y')

        save_filename = ['PP_PnP_',relevant_condition,'_',rel_save_name,'_',opacity,'_',modality_name,'_confid_score_',num2str(confidence_score)];

    else
        
        save_filename = ['PP_PnP_',relevant_condition,'_',rel_save_name,'_',opacity,'_',modality_name,'_confid_score_',num2str(confidence_score),'_no_blink_rejections'];

    end
    
end

%Define and make save directory
save_dir = fullfile('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Group Analysis\EyeLink Analysis\Timecourses with Stats',save_filename);
mkdir(save_dir)

%Figure name 
if isequal(trial_type,'pt')
    
    figure_name = 'PP_vs_PnP';
   
elseif isequal(trial_type,'ct')
    
    figure_name = 'CP_vs_CnP';
    
end

%% Load Group EyeLink Data 

%Load group Eyelink data
cd(eyelink_dir)
load(eyelink_filename)

%Rename data to generic variable 
if isequal(trial_type,'pt')

   group_P_pupil_data = group_PP_pupil_data;
   group_nP_pupil_data = group_PnP_pupil_data;
   
   group_P_blink_data = group_PP_blink_data;
   group_nP_blink_data = group_PnP_blink_data;
   
   group_P_microsac_data = group_PP_microsac_data;
   group_nP_microsac_data = group_PnP_microsac_data;
    
elseif isequal(trial_type,'ct')
    
   group_P_pupil_data = group_CP_pupil_data;
   group_nP_pupil_data = group_CnP_pupil_data;
   
   group_P_blink_data = group_CP_blink_data;
   group_nP_blink_data = group_CnP_blink_data;
   
   group_P_microsac_data = group_CP_microsac_data;
   group_nP_microsac_data = group_CnP_microsac_data;
    
end

%% Preprocess Data for Visualization

%Define stim start time
stim_time = 6001;

%Baseline Pupil Data

%Define baseline period
baseline_window = 5001:6000;

%Create a baseline matrix across subjects
P_baseline_matrix = nanmean(group_P_pupil_data(baseline_window,:),1)';
nP_baseline_matrix = nanmean(group_nP_pupil_data(baseline_window,:),1)';

%Subtract baseline matrix from group pupil data
P_pupil_baselined = group_P_pupil_data - shiftdim(repmat(P_baseline_matrix,1,12000),1);
nP_pupil_baselined = group_nP_pupil_data - shiftdim(repmat(nP_baseline_matrix,1,12000),1);

%Bin Microsaccade Data

%Moving average bin size
bin_size = 500;

%Moving average over microsaccade time points
P_microsaccade_smoothed = movmean(group_P_microsac_data,bin_size);
nP_microsaccade_smoothed = movmean(group_nP_microsac_data,bin_size);

%% Define Plot Color

if isequal(relevant_condition, 'relevant')
    
    P_color = 'b';
    nP_color = 'r';
    
elseif isequal(relevant_condition, 'irrelevant') || isequal(relevant_condition, 'rel_irrel')
    
    P_color = 'g';
    nP_color = 'y';
    
end

%% Pupil Diameter

%Load stats
cd(stat_dir)
load(['pupil_timecourse_cluster_',figure_name,'.mat']);

%Define sig line y-axis
current_yaxis = 0.14;
%current_yaxis = 0.16;

%Setup figure
pupil_plot = figure
hold on

title(['Pupil - Right Eye - ', relevant_location, ' - ', modality_name],'Interpreter','none')
xlabel('Time (ms)')
ylabel('Baselined Diameter (mm)')

%Define dimensions
ylim([-0.15 0.15])
%ylim([-0.125 0.175])
%xlim([-1000 5000])
xlim([-4000 5000])

%Tick marks 
yticks([-0.15:0.05:0.15])
%yticks([-0.125:0.025:0.175])

%Timevector/time period (Note: the timevector includes 0 as the stim time)
%timevector = [-1000:4999];
%time_period = [5001:11000];
timevector = [-5000:4999];
time_period = [1001:11000];

%Plot stim onset line
plot([0 0],[-0.15,0.15],'k')

%Plot P and nP mean baselined timecourses
P_plot = plot(timevector, nanmean(P_pupil_baselined(time_period,:),2),P_color)
nP_plot = plot(timevector, nanmean(nP_pupil_baselined(time_period,:),2),nP_color)

%Plot CP SEM timecoureses
plot(timevector, nanmean(P_pupil_baselined(time_period,:),2) + ...
    (std(P_pupil_baselined(time_period,:),0,2)/sqrt(size(P_pupil_baselined,2))),['--',P_color])

plot(timevector, nanmean(P_pupil_baselined(time_period,:),2) - ...
    (std(P_pupil_baselined(time_period,:),0,2)/sqrt(size(P_pupil_baselined,2))),['--',P_color])

%Plot CnP SEM timecoureses
plot(timevector, nanmean(nP_pupil_baselined(time_period,:),2) + ...
    (std(nP_pupil_baselined(time_period,:),0,2)/sqrt(size(nP_pupil_baselined,2))),['--',nP_color])

plot(timevector, nanmean(nP_pupil_baselined(time_period,:),2) - ...
    (std(nP_pupil_baselined(time_period,:),0,2)/sqrt(size(nP_pupil_baselined,2))),['--',nP_color])

%Plot statistically significant times

%If there are significant times points
if not(isempty(sig_time_pts))

    %Find significant break points
    sig_breaks = find(diff(sig_time_pts) > 1);

    %If there breaks in sig times
    if not(isempty(sig_breaks))
        
        %Loop over discontinuous points
        for breaks = 1:length(sig_breaks)

            %Current break point
            break_point = sig_breaks(breaks);

            %First break
            if breaks == 1 

                plot((sig_time_pts(1:break_point)-stim_time),ones(1,length(sig_time_pts(1:break_point)))*current_yaxis, 'Color', 'k', 'LineWidth',2)

                %If only one break point also plot to the end of
                %sig times
                if length(sig_breaks) == 1

                    %Between last break and end
                    plot([sig_time_pts(break_point+1:length(sig_time_pts))]-stim_time,ones(1,length(sig_time_pts(break_point+1:length(sig_time_pts))))*current_yaxis, 'Color', 'k', 'LineWidth',2)

                end

            %Last break
            elseif breaks == length(sig_breaks)

                %Between penultimate and final break
                plot([sig_time_pts(sig_breaks(breaks-1)+1):sig_time_pts(break_point)]-stim_time,...
                    ones(1,length((sig_time_pts(sig_breaks(breaks-1)+1):sig_time_pts(break_point))))*current_yaxis, 'Color', 'k', 'LineWidth',2)

                %Between last break and end
                plot([sig_time_pts(break_point+1:length(sig_time_pts))]-stim_time,ones(1,length(sig_time_pts(break_point+1:length(sig_time_pts))))*current_yaxis, 'Color', 'k', 'LineWidth',2)

            %Middle breaks
            else

                plot([sig_time_pts(sig_breaks(breaks-1)+1):sig_time_pts(break_point)]-stim_time,...
                    ones(1,length((sig_time_pts(sig_breaks(breaks-1)+1):sig_time_pts(break_point))))*current_yaxis, 'Color', 'k', 'LineWidth',2)

            end

        end
    
    %If the sig time points are continuous
    else

      %Plot sig line
      plot(sig_time_pts-stim_time, ones(1,length(sig_time_pts))*current_yaxis, 'Color', 'k', 'LineWidth',2)
    
    end
    
end

%Plot figure legend
legend([P_plot,nP_plot],['Perceived (n = ', num2str(size(P_pupil_baselined,2)),')'], ...
    ['Not Perceived (n = ', num2str(size(nP_pupil_baselined,2)),')'],'Location','best')

%Save figure
cd(save_dir)
savefig(pupil_plot,['Pupil_timecourse_',figure_name,'.fig'])
%saveas(pupil_plot,['Pupil_timecourse_',figure_name,'.eps'],'eps')

close

%% Blink Rate 

%Load stats
cd(stat_dir)
load(['blink_timecourse_cluster_',figure_name,'.mat']);

%Define sig line y-axis
current_yaxis = 0.21;
%current_yaxis = 0.24;

%Setup figure
blink_plot = figure
hold on

title(['Blink - Right Eye - ', relevant_location, ' - ', modality_name],'Interpreter','none')
xlabel('Time (ms)')
ylabel('Blink occurrence')

%Define dimensions
ylim([0 0.22])
%ylim([0 0.27])
xlim([-1000 5000])
%xlim([-4000 5000])

yticks([0.02:0.04:0.22])
%yticks([0.02:0.04:0.26]) %Used for non-combined datasets 

%Timevector/time period (Note: the timevector includes 0 as the stim time)
timevector = [-1000:4999];
time_period = [5001:11000];
%timevector = [-5000:4999];
%time_period = [1001:11000];

%Plot stim onset line
plot([0 0],[0,0.22],'k')

%Plot CP and CnP mean baselined timecourses (Note: multiple by 100 if you
%want percent, otherwise the value represents rate)
P_plot = plot(timevector, nanmean(group_P_blink_data(time_period,:),2),P_color)
nP_plot = plot(timevector, nanmean(group_nP_blink_data(time_period,:),2),nP_color)

%Plot CP SEM timecoureses
plot(timevector, nanmean(group_P_blink_data(time_period,:),2) + ...
    (std(group_P_blink_data(time_period,:),0,2)/sqrt(size(group_P_blink_data,2))),['--',P_color])

plot(timevector, nanmean(group_P_blink_data(time_period,:),2) - ...
    (std(group_P_blink_data(time_period,:),0,2)/sqrt(size(group_P_blink_data,2))),['--',P_color])

%Plot CnP SEM timecoureses
plot(timevector, nanmean(group_nP_blink_data(time_period,:),2) + ...
    (std(group_nP_blink_data(time_period,:),0,2)/sqrt(size(group_nP_blink_data,2))),['--',nP_color])

plot(timevector, nanmean(group_nP_blink_data(time_period,:),2) - ...
    (std(group_nP_blink_data(time_period,:),0,2)/sqrt(size(group_nP_blink_data,2))),['--',nP_color])

%Plot statistically significant times

%If there are significant times points
if not(isempty(sig_time_pts))

    %Find significant break points
    sig_breaks = find(diff(sig_time_pts) > 1);

    %If there breaks in sig times
    if not(isempty(sig_breaks))
        
        %Loop over discontinuous points
        for breaks = 1:length(sig_breaks)

            %Current break point
            break_point = sig_breaks(breaks);

            %First break
            if breaks == 1 

                plot((sig_time_pts(1:break_point)-stim_time),ones(1,length(sig_time_pts(1:break_point)))*current_yaxis, 'Color', 'k', 'LineWidth',2)

                %If only one break point also plot to the end of
                %sig times
                if length(sig_breaks) == 1

                    %Between last break and end
                    plot([sig_time_pts(break_point+1:length(sig_time_pts))]-stim_time,ones(1,length(sig_time_pts(break_point+1:length(sig_time_pts))))*current_yaxis, 'Color', 'k', 'LineWidth',2)

                end

            %Last break
            elseif breaks == length(sig_breaks)

                %Between penultimate and final break
                plot([sig_time_pts(sig_breaks(breaks-1)+1):sig_time_pts(break_point)]-stim_time,...
                    ones(1,length((sig_time_pts(sig_breaks(breaks-1)+1):sig_time_pts(break_point))))*current_yaxis, 'Color', 'k', 'LineWidth',2)

                %Between last break and end
                plot([sig_time_pts(break_point+1:length(sig_time_pts))]-stim_time,ones(1,length(sig_time_pts(break_point+1:length(sig_time_pts))))*current_yaxis, 'Color', 'k', 'LineWidth',2)

            %Middle breaks
            else

                plot([sig_time_pts(sig_breaks(breaks-1)+1):sig_time_pts(break_point)]-stim_time,...
                    ones(1,length((sig_time_pts(sig_breaks(breaks-1)+1):sig_time_pts(break_point))))*current_yaxis, 'Color', 'k', 'LineWidth',2)

            end

        end
    
    %If the sig time points are continuous
    else

      %Plot sig line
      plot(sig_time_pts-stim_time, ones(1,length(sig_time_pts))*current_yaxis, 'Color', 'k', 'LineWidth',2)
    
    end
    
end

%Save figure
cd(save_dir)
savefig(blink_plot,['Blink_timecourse_',figure_name,'.fig'])
%saveas(blink_plot,['Blink_timecourse_',figure_name,'.eps'],'eps')

close

%% Microsaccade Rate

%Load stats
cd(stat_dir)
load(['microsac_timecourse_cluster_',figure_name,'.mat']);

%Define sig line y-axis
current_yaxis = 0.017;

if isequal(modality_name,'EEG')

    current_yaxis = 0.027;
  
elseif isequal(modality_name,'MRI')
    
    current_yaxis = 0.027;
end

%Setup figure
microsaccade_plot = figure
hold on

title(['Microsaccade - Right Eye - ', relevant_location, ' - ', modality_name],'Interpreter','none')
xlabel('Time (ms)')
ylabel('Microsaccade occurence')

%Define dimensions
ylim([0.005 0.018])
%xlim([-1000 5000])
xlim([-4000 5000])

%yticks([0.006:0.002:0.018])

if isequal(modality_name,'EEG')
    
    %yticks([0.006:0.002:0.018])
    %ylim([0.006 0.018])
    yticks([0.002:0.002:0.028])
    ylim([0.002 0.028])
    
elseif isequal(modality_name,'MRI')
    
    yticks([0.002:0.002:0.028])
    ylim([0.002 0.028])

end

%Timevector/time period (Note: the timevector includes 0 as the stim time)
timevector = [-1000:4999];
time_period = [5001:11000];
%timevector = [-5000:4999];
%time_period = [1001:11000];

%Plot stim onset line
plot([0 0],[0,0.018],'k')

%Plot CP and CnP mean baselined timecourses
P_plot = plot(timevector, nanmean(P_microsaccade_smoothed(time_period,:),2),P_color);
nP_plot = plot(timevector, nanmean(nP_microsaccade_smoothed(time_period,:),2),nP_color)

%Plot CP SEM timecoureses
plot(timevector, nanmean(P_microsaccade_smoothed(time_period,:),2) + ...
    (std(P_microsaccade_smoothed(time_period,:),0,2)/sqrt(size(P_microsaccade_smoothed,2))),['--',P_color])

plot(timevector, nanmean(P_microsaccade_smoothed(time_period,:),2) - ...
    (std(P_microsaccade_smoothed(time_period,:),0,2)/sqrt(size(P_microsaccade_smoothed,2))),['--',P_color])

%Plot CnP SEM timecoureses
plot(timevector, nanmean(nP_microsaccade_smoothed(time_period,:),2) + ...
    (std(nP_microsaccade_smoothed(time_period,:),0,2)/sqrt(size(nP_microsaccade_smoothed,2))),['--',nP_color])

plot(timevector, nanmean(nP_microsaccade_smoothed(time_period,:),2) - ...
    (std(nP_microsaccade_smoothed(time_period,:),0,2)/sqrt(size(nP_microsaccade_smoothed,2))),['--',nP_color])

%Plot statistically significant times

%If there are significant times points
if not(isempty(sig_time_pts))

    %Find significant break points
    sig_breaks = find(diff(sig_time_pts) > 1);

    %If there breaks in sig times
    if not(isempty(sig_breaks))
        
        %Loop over discontinuous points
        for breaks = 1:length(sig_breaks)

            %Current break point
            break_point = sig_breaks(breaks);

            %First break
            if breaks == 1 

                plot((sig_time_pts(1:break_point)-stim_time),ones(1,length(sig_time_pts(1:break_point)))*current_yaxis, 'Color', 'k', 'LineWidth',2)

                %If only one break point also plot to the end of
                %sig times
                if length(sig_breaks) == 1

                    %Between last break and end
                    plot([sig_time_pts(break_point+1:length(sig_time_pts))]-stim_time,ones(1,length(sig_time_pts(break_point+1:length(sig_time_pts))))*current_yaxis, 'Color', 'k', 'LineWidth',2)

                end

            %Last break
            elseif breaks == length(sig_breaks)

                %Between penultimate and final break
                plot([sig_time_pts(sig_breaks(breaks-1)+1):sig_time_pts(break_point)]-stim_time,...
                    ones(1,length((sig_time_pts(sig_breaks(breaks-1)+1):sig_time_pts(break_point))))*current_yaxis, 'Color', 'k', 'LineWidth',2)

                %Between last break and end
                plot([sig_time_pts(break_point+1:length(sig_time_pts))]-stim_time,ones(1,length(sig_time_pts(break_point+1:length(sig_time_pts))))*current_yaxis, 'Color', 'k', 'LineWidth',2)

            %Middle breaks
            else

                plot([sig_time_pts(sig_breaks(breaks-1)+1):sig_time_pts(break_point)]-stim_time,...
                    ones(1,length((sig_time_pts(sig_breaks(breaks-1)+1):sig_time_pts(break_point))))*current_yaxis, 'Color', 'k', 'LineWidth',2)

            end

        end
    
    %If the sig time points are continuous
    else

      %Plot sig line
      plot(sig_time_pts-stim_time, ones(1,length(sig_time_pts))*current_yaxis, 'Color', 'k', 'LineWidth',2)
    
    end
    
end

%Save figure
cd(save_dir)
savefig(microsaccade_plot,['Microsaccade_timecourse_',figure_name,'.fig'])
%saveas(microsaccade_plot,['Microsaccade_timecourse_',figure_name,'.eps'],'eps')

close
