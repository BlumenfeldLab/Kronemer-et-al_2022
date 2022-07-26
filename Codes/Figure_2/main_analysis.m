%% Kronemer et al., 2021 
% This script preprocess low-density EEG and icEEG data and 
% plots the ldEEG v. icEEG ERP trace for Figure 3. 
% 
% Written by Thomas Xin 
% Sep. 4, 2021 
%
%Written by: Thomas Xin  
%Date: 9/4/2021

%% add helper functions to path
addpath('R:\RNS Study\Analysis\Analysis Code\Temporal Analysis'); 
addpath('R:\RNS Study\Analysis\Analysis Code\Temporal Analysis\helper functions');
%% run the script to read RNS locations
% only need to run once 

% run on server 
% process_451NH ; 

% run locally
% read_RNS_location ;

%% Set parameters

%Sampling rate
sampling_rate = 250; 

%Subject IDs
subjects = {'471MH','489RD','491GS','490KB','535BP','536JD'};

% all eeg names in low-density EEG
reference_map = lower({'Fp1','Fpz','Fp2','F7','F3','Fz','F4','F8','C3','Cz','C4', 'P3','Pz','P4','O1','O2'});

% main channels of analysis 
main_ch = {'oz','cz','pz','fz'};
main_ch_idx = [17,10,13,6]; 

%Reference electrode
reference_RNS_channel = {'CM','MB','CM&MB','noCMMB'}; 

% polarity by visual inspection. 1 - positive; 0 - negative; NaN - below
% threshold - this would be consistent with Fig. 2(D) 
% if not sure about polarity. use NaN as placeholder. Run the script and
% add polarity index later. 
polarity_idx = [1,0,1,1,0,0,1,1,0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];

set(0, 'DefaultLineLineWidth', 1.3);
set(0,'defaultAxesFontSize',10); 

%% Load EEG and RNS data 

RNS_EEG_Manager = DataManager(sampling_rate,subjects,reference_map);
RNS_EEG_Manager = RNS_EEG_Manager.loadData(); 
RNS_EEG_Manager = RNS_EEG_Manager.addOz(); % extrapolate Oz channel 

%% Load 451NH 

RNS_451NH_CP = load('\\server2.med.yale.internal\Data7\RNS Study\Analysis\Subject Analysis\451NH\icEEG Analysis\cp_CM_channel_epoch_CP.mat','cp_data_plot'); 
RNS_451NH_CP = RNS_451NH_CP.cp_data_plot(1:2,:,:); % only load the first two channels of 451 

RNS_451NH_CnP = load('\\server2.med.yale.internal\Data7\RNS Study\Analysis\Subject Analysis\451NH\icEEG Analysis\CM_channel_epoch_CnP.mat','cnp_data_plot'); 
RNS_451NH_CnP = RNS_451NH_CnP.cnp_data_plot(1:2,:,:); % only load the first two channels of 451 

%% Load RNS data 

RNS_CH_MAP = '\\server2.med.yale.internal\Data7\RNS Study\Analysis\Group Analysis\RNS_channel_table.mat' ;
RNS_CP = RNSAnalyzer(RNS_EEG_Manager.RNS_CP,RNS_CH_MAP,sampling_rate,subjects,reference_RNS_channel,RNS_451NH_CP,polarity_idx); 
RNS_CnP = RNSAnalyzer(RNS_EEG_Manager.RNS_CnP,RNS_CH_MAP,sampling_rate,subjects,reference_RNS_channel,RNS_451NH_CnP,ones(size(polarity_idx)));

%% Fig.2 Panel D - inspect polarity and ampltiude 
% visualize all subject-channel pairs ranked by post-stimulus absolute
% amplitude

RNS_CP = RNS_CP.plotSNRBySubject([-25,25],'amplitude');

%% Flip channels with negative first peak in ERP 

RNS_EEG_Manager.RNS_CP = RNS_CP.flipSourceData(); 
RNS_EEG_Manager.RNS_CnP  = RNS_CnP.flipSourceData(); 

%% Lowpass filter RNS and EEG data 

RNS_EEG_Manager = RNS_EEG_Manager.filterLowPass(); 
RNSEEG_CP = EEGRNSAnalyzer( RNS_EEG_Manager.EEG_CP, RNS_EEG_Manager.RNS_CP,RNS_CP,RNS_EEG_Manager.EEG_reference_map,reference_RNS_channel,RNS_EEG_Manager);
RNSEEG_CnP = EEGRNSAnalyzer( RNS_EEG_Manager.EEG_CnP,RNS_EEG_Manager.RNS_CnP,RNS_CnP,RNS_EEG_Manager.EEG_reference_map,reference_RNS_channel,RNS_EEG_Manager);

%% Plot all erps 

RNSEEG_CP = RNSEEG_CP.plotEEGRNSTrace('Lowpass Filtered CP','Voltage (uv)',[-6,8],'ERP') ;
RNSEEG_CnP = RNSEEG_CnP.aggregateEEGRNSTrace('Lowpass Filtered CnP','Voltage (uv)',[-6,8],'ERP') ;
RNSEEG_CP.plotEEGRNSTraceCombined(RNSEEG_CnP,1,'ERP Trace CP&CnP','Voltage (uv)',[-6,10],'ERP','both',main_ch,false) ;

%% Plot arrival latency  

RNSEEG_CP = RNSEEG_CP.getArrivalLatency(RNS_EEG_Manager.epoch_length,main_ch);

%% Save for latency analysis 
save('R:\RNS Study\Analysis\Group Analysis\Correlation Analysis\rns_cp_small.mat','RNSEEG_CP');

%% Output for cluster perm and bootstrap 

group_rns_cp = RNSEEG_CP.RNS_data; 
group_rns_cnp = RNSEEG_CnP.RNS_data; 
group_eeg_cp = RNSEEG_CP.EEG_data;  
group_eeg_cnp = RNSEEG_CnP.EEG_data; 

save('group_rns.mat','group_rns_cp', 'group_rns_cnp');
save('group_eeg.mat','group_eeg_cp', 'group_eeg_cnp');
