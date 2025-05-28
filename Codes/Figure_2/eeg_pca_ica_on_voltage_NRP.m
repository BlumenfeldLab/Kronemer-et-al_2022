%% Run PCA and ICA

%Written by: Sharif I. Kronemer and Mark Aksen
%Date: 7/23/2019
%Modified: 1/19/2021

function [ICA_epochs_components_removed, remove_components] = eeg_pca_ica_on_voltage_NRP(ID, event_type, ICA_epochs, Photo_dir, preprocessed_dir, sampling_rate, num_components, eeglab_prepared_blank)

%Rename EEGLab Structure
ICA_eeglab_prepared_blank = eeglab_prepared_blank;

%Setup EEGLab Structure
ICA_eeglab_prepared_blank.data = ICA_epochs; %Define data
ICA_eeglab_prepared_blank.nbchan = size(ICA_epochs,1); %Define number of channels
ICA_eeglab_prepared_blank.xmax = size(ICA_epochs,2)/sampling_rate; %Define end time
ICA_eeglab_prepared_blank.times = 0:size(ICA_epochs,2)-1; %Define time points
ICA_eeglab_prepared_blank.pnts = size(ICA_epochs,2); %Define number of samples
ICA_eeglab_prepared_blank.trials = size(ICA_epochs,3); %Define number of trials
ICA_eeglab_prepared_blank.ref = 'averef'; %Define reference type

% %Load channel location data
% cd(Photo_dir)
% chanlocs = readlocs([ID, '_GPS.sfp']);
% 
% %Remove nasion points
% chanlocs([1:3,260]) = []; 
% ICA_eeglab_prepared_blank.chanlocs = chanlocs; %Update channel location

%Run PCA - ICA    
[weightsrecut,~] = pop_runica(ICA_eeglab_prepared_blank,'icatype','runica','pca',num_components);

%Assign proper channel locations
weightsrecut.chanlocs = ICA_eeglab_prepared_blank.chanlocs; %chanlocs;

%Display ICA components
display = pop_selectcomps(weightsrecut);

%Select which components to remove (Need to type into command line)
prompt_ICA = 'Components you want to remove (e.g., [1, 2, 4]; or "none"): ';
remove_components = input(prompt_ICA);

%Save ICA component topoplots and matrix of removed components 
cd(preprocessed_dir)

if isequal(event_type,'faces')
    
    savefig('rel_irrel_faces_epochs_voltage_ICA_components_topoplot.fig')

elseif isequal(event_type, 'questions_buttons')
    
    savefig('questions_buttons_epochs_voltage_ICA_components_topoplot.fig')
    
elseif isequal(event_type, 'faces_post_jitter')
    
    savefig('faces_post_jitter_epochs_voltage_ICA_components_topoplot.fig')

end

%Close the windows that are open
close all

    %Remove ICA components
    if isequal(remove_components, 'none')

        disp('***Not removing ICA components***');

        %Rename variable for output of function
        ICA_epochs_components_removed = ICA_epochs;

    else

        disp('***Removing ICA components***');

        %Remove components from data (inputs: data, matrix of components, display channel pre-post ICA removal, remove components)
        weightsrecut_subcomp = pop_subcomp(weightsrecut,remove_components,1,0);

        %ICA component removed epochs
        ICA_epochs_components_removed = weightsrecut_subcomp.data;

    end

end