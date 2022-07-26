%% Plot Pupil, Blink, and Microsaccade Timecourses - PP/PnP Over Multiple Confidence Scores

%This code will load the group Eyelink data and the group eyelink data
%permutation statistical test results and plot the eyelink timecourses 

%Written by: Sharif I. Kronemer
%Date: 7/5/2021
%Last Modified: 5/15/2022

clear
close all

%% Prompt

%Relevant condition
prompt_1 = 'Relevant condition [relevant, irrelevant, rel_irrel]: ';
relevant_condition = input(prompt_1, 's');

%Location relevant
prompt_2 = 'No-report relevant location [c, q, combine]: ';
relevant_location = input(prompt_2, 's');

%Opacity
prompt_3 = 'Stimulus opacity [threshold, blank]: ';
opacity = input(prompt_3, 's');

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

%Confidence Score
confidence_score_matrix = [0:0.25:2];

%Irrelevant condition
if isequal(relevant_condition, 'irrelevant')
    
    %Center Relevant/Irrelevant
    if isequal(relevant_location, 'c')

        relevant_location = 'Center Irrelevant';
        rel_save_name = 'cent';

    %Quadrant Relevant
    elseif isequal(relevant_location, 'q')

        relevant_location = 'Quadrant Irrelevant';
        rel_save_name = 'quad';

    %Combine Center and Quadrant Relevant   
    elseif isequal(relevant_location, 'combine')

        relevant_location = 'Center and Quadrant Irrelevant';      
        rel_save_name = 'cent_quad';

    end

%Relevant condition
elseif isequal(relevant_condition, 'relevant')
    
    %Center Relevant/Irrelevant
    if isequal(relevant_location, 'c')

        relevant_location = 'Center Relevant';
        rel_save_name = 'cent';

    %Quadrant Relevant
    elseif isequal(relevant_location, 'q')

        relevant_location = 'Quadrant Relevant';
        rel_save_name = 'quad';

    %Combine Center and Quadrant Relevant   
    elseif isequal(relevant_location, 'combine')

        relevant_location = 'Center and Quadrant Relevant';      
        rel_save_name = 'cent_quad';

    end
    
%Relevant and Irrelevant
elseif isequal(relevant_condition, 'rel_irrel')
    
       %Center Relevant/Irrelevant
    if isequal(relevant_location, 'c')

        relevant_location = 'Center Relevant and Irrelevant';
        rel_save_name = 'cent';

    %Quadrant Relevant
    elseif isequal(relevant_location, 'q')

        relevant_location = 'Quadrant Relevant and Irrelevant';
        rel_save_name = 'quad';

    %Combine Center and Quadrant Relevant   
    elseif isequal(relevant_location, 'combine')

        relevant_location = 'Center and Quadrant Relevant and Irrelevant';      
        rel_save_name = 'cent_quad';

    end
    
end

%Figure name 
figure_name = 'PP_vs_PnP';

%Eyelink directory
eyelink_dir = 'S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Group Data\EyeLink';

%Save Filename
save_filename = ['PP_PnP_',relevant_condition,'_',rel_save_name,'_',opacity,'_',modality_name,'_all_confid_score'];

%Save dir
save_dir = fullfile('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Group Analysis\EyeLink Analysis\Timecourses with Stats',save_filename);
mkdir(save_dir)

%% Setup Figures

set(0, 'DefaultFigureRenderer', 'painters');
 
%PUPIL FIGURE

%Setup figure
pupil_plot = figure(1)
hold on

title(['Pupil - Right Eye - ', relevant_location, ' - ', modality_name,' - Scores 0-2'],'Interpreter','none')
xlabel('Time (ms)')
ylabel('Baselined Diameter (mm)')

%Define dimensions
ylim([-0.45 0.45])
xlim([-1000 5000])

%Tick marks 
%yticks([-0.125:0.025:0.175])

%Timevector/time period (Note: the timevector includes 0 as the stim time)
timevector = [-1000:4999];
time_period = [5001:11000];

%Plot stim onset line
plot([0 0],[-0.45,0.45],'k')

%BLINK FIGURE

%Setup figure
blink_plot = figure(2)
hold on

title(['Blink - Right Eye - ', relevant_location, ' - ', modality_name,' - Scores 0-2'],'Interpreter','none')
xlabel('Time (ms)')
ylabel('Blink occurrence')

%Define dimensions
ylim([0 0.65])
xlim([-1000 5000])

%yticks([0.02:0.04:0.26]) %Used for non-combined datasets 

%Timevector/time period (Note: the timevector includes 0 as the stim time)
timevector = [-1000:4999];
time_period = [5001:11000];

%Plot stim onset line
plot([0 0],[0,0.65],'k')

%MICROSACCADE FIGURE

%Setup figure
microsaccade_plot = figure(3)
hold on

title(['Microsaccade - Right Eye - ', relevant_location, ' - ', modality_name,' - Scores 0-2'],'Interpreter','none')
xlabel('Time (ms)')
ylabel('Microsaccade occurence')

%Define dimensions
xlim([-1000 5000])
ylim([0 0.055])
yticks([0:0.005:0.055])

% if isequal(modality_name,'EEG')
%     
%     yticks([0.002:0.002:0.028])
%     ylim([0.002 0.028])
%     
% elseif isequal(modality_name,'MRI')
%     
%     yticks([0.002:0.002:0.028])
%     ylim([0.002 0.028])
% 
% end

%Timevector/time period (Note: the timevector includes 0 as the stim time)
timevector = [-1000:4999];
time_period = [5001:11000];

%Plot stim onset line
plot([0 0],[0,0.055],'k')

%Color vector
color_vector = {'r','b','g','y','c','m','k',[0.5,0.18,0.55],[0.85,0.3,0.1]};
  
%Loop over confidence score
for score = 1:length(confidence_score_matrix)
    
    %Current confidence score
    confidence_score = confidence_score_matrix(score);
    
    disp(['Running confidence score ', num2str(confidence_score)])

    %Current color
    current_color = color_vector{score};

    P_color = current_color;
    nP_color = current_color;
    
    %% Directories 
    
    %Eyelink filename
    eyelink_filename = ['Group_eyelink_',relevant_condition,'_PP_PnP_',rel_save_name,'_',opacity,'_',modality_name,'_score_thres_',num2str(confidence_score),'_data.mat'];

    %Statistical result directory
    stat_dir = fullfile('S:\HNCT_NRP_Study\HNCT No Report Paradigm\Analysis\Group Analysis\EyeLink Analysis\Permutation Analysis',...
        ['PP_PnP_',relevant_condition,'_',rel_save_name,'_',opacity,'_',modality_name,'_confid_score_',num2str(confidence_score)]);

    %% Load Group EyeLink Data 

    %Load group Eyelink data
    cd(eyelink_dir)
    load(eyelink_filename)
    
    group_P_pupil_data = group_PP_pupil_data;
    group_nP_pupil_data = group_PnP_pupil_data;

    group_P_blink_data = group_PP_blink_data;
    group_nP_blink_data = group_PnP_blink_data;

    group_P_microsac_data = group_PP_microsac_data;
    group_nP_microsac_data = group_PnP_microsac_data;

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

    %% Pupil Diameter
    
    %pupil_plot
    figure(1)

    %Load stats
    cd(stat_dir)
    load(['pupil_timecourse_cluster_',figure_name,'.mat']);

    %Define sig line y-axis
    %current_yaxis = 0.14;
    current_yaxis = 0.16;

    %Plot P and nP mean baselined timecourses
    P_plot = plot(timevector, nanmean(P_pupil_baselined(time_period,:),2),'color',P_color)
    nP_plot = plot(timevector, nanmean(nP_pupil_baselined(time_period,:),2),'color',nP_color)

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

    %% Blink Rate 

    figure(2)
    
    %Load stats
    cd(stat_dir)
    load(['blink_timecourse_cluster_',figure_name,'.mat']);

    %Define sig line y-axis
    %current_yaxis = 0.21;
    current_yaxis = 0.24;

    %Plot CP and CnP mean baselined timecourses
    P_plot = plot(timevector, nanmean(group_P_blink_data(time_period,:),2),'color',P_color)
    nP_plot = plot(timevector, nanmean(group_nP_blink_data(time_period,:),2),'color',nP_color)

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

    %% Microsaccade Rate

    figure(3)
    
    %Load stats
    cd(stat_dir)
    load(['microsac_timecourse_cluster_',figure_name,'.mat']);

    %Define sig line y-axis
    %current_yaxis = 0.017;

    if isequal(modality_name,'EEG')

        current_yaxis = 0.027;

    elseif isequal(modality_name,'MRI')

        current_yaxis = 0.027;
    end

    %Plot CP and CnP mean baselined timecourses
    P_plot = plot(timevector, nanmean(P_microsaccade_smoothed(time_period,:),2),'color',P_color);
    nP_plot = plot(timevector, nanmean(nP_microsaccade_smoothed(time_period,:),2),'color',nP_color)

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

end

%Save figure
cd(save_dir)
savefig(pupil_plot,['Pupil_timecourse_',figure_name,'.fig'])

%Save figure
cd(save_dir)
savefig(blink_plot,['Blink_timecourse_',figure_name,'.fig'])

%Save figure
cd(save_dir)
savefig(microsaccade_plot,['Microsaccade_timecourse_',figure_name,'.fig'])
