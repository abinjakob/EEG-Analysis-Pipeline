% % SSVEP EEG DATA PROCESSING MAIN SCRIPT 
% ------------------------------------
% The code performs the pre-processing and processing pipeline for 
% the artefact corrected EEG data 
% 
% Pre-processing steps:
% 1. Filtering (Low-pass and High-pass filters)
% 2. Rereference to CAR
% 3. Epoching data
% 4. Baseline Correction
% 5. Reject artefactual epochs
%
% Data-procesing steps:
% 6. Sub-epoching the data 
% 7. Calculating PSD for sub-epoched data
% 8. Plot PSD for selected channel and condition
% 
% the pre-processed data is then stored to the given folder
% 
% Author: Abin Jacob
% Date  : 28/11/2023

%% start fresh 
clear; clc; close all;

%% define project folders

% path to the main folder 
MAINPATH = '/Users/abinjacob/Documents/02. NeuroCFN/Research Module/Paradigm/Analysis Codes/EEG Data Analysis';
% path for the rawdata
DATAPATH = fullfile(MAINPATH, 'ica_corrected_data','ICA_cleaned',filesep);
% path to save output
PATHOUT = fullfile(MAINPATH, 'pre_processed_data',filesep);

% add EEGLab to matlab path
addpath('L:\Cloud\SW\eeglab2023.1');
% open EEGLab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% set current directory 
cd(MAINPATH);

%% SET PARAMETERS FOR DATA PROCESSING

% load ICA-clean file
%EEG = pop_loadset('/Users/abinjacob/Documents/02. NeuroCFN/Research Module/EEG Analysis Scripts/temp_rawdata/SSVEP_pilot2_ICAcleaned.set');

% filtering 
% high-pass filter 
HP = 0.1;                       % cut-off
HP_order = 826;                 % filter order    
% low-pass filter  
LP = 50;                        % cut-off
LP_order = 776;                 % filter order 

% epoching
% event markers 
events = {'stim_L20','stim_L15','stim_R20','stim_R15'};
epoch_start = -0.2;            
epoch_end = 4;               
% baseline correction
% defining baseline for baseline correcion
baseline = [epoch_start*EEG.srate 0];   

% reject artefactual epochs 
PRUNE = 4;

%% pre-processing 

% low-pass filter
EEG = pop_firws(EEG, 'fcutoff', LP, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', LP_order);
% high-pass filter
EEG = pop_firws(EEG, 'fcutoff', HP, 'ftype', 'highpass', 'wtype', 'hamming', 'forder', HP_order);
% re-referencing to CAR
EEG = pop_reref(EEG, [], 'refstate',0);

% removing unnecessary event marker
event_pos = 1;      % position counter for the events other than stim onset
event_idx = [];     % array to store the index of the event other than stim onset
% loop over events 
for idx = 1: length(EEG.event)
    if ~ strcmp(EEG.event(idx).type, events)
        event_idx(event_pos) = idx;
        event_pos = event_pos +1;
    end
end 
% remove events which are not stim onset from the data
EEG = pop_editeventvals(EEG, 'delete', event_idx);
EEG = eeg_checkset(EEG);

% epoching 
EEG = pop_epoch(EEG, events, [epoch_start epoch_end], 'newname', 'MI_pilot_epoched','epochinfo', 'yes');
% reject artefactual epochs 
% joint probability-based artifact rejection (joint prob. > PRUNE (SD))
EE_temp = pop_jointprob(EEG, 1, [1:EEG.nbchan], PRUNE, PRUNE, 0, 1, 0);
EEG = eeg_checkset(EEG);
% baseline correction
EEG = pop_rmbase(EEG, baseline);
EEG = eeg_checkset(EEG);

%% epcohing and calculating ERD for each condition

% loop over events 
for iEvent = 1:length(events)
    % epoching  
    EEG_temp = pop_selectevent(EEG, 'type', events{iEvent},'renametype', events{iEvent}, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
    
    % window size for pwelch 
    window_length = EEG_temp.pnts;

    % calculating psd and storing
    [event(iEvent).psd, event(iEvent).f] = calc_psd(EEG_temp, window_length);
    % setting event type
    event(iEvent).eventtype = events{iEvent};
end 

%% plotting PSD

% params for plotting
event_sel = 1;       % event to plot
channel_sel = 18;    % channel to plot
stims = {'Left 20Hz', 'Left 15Hz', 'Right 20Hz', 'Right 15Hz'};

% plotting PSD for selected channel
figure;
% plot(event(event_sel).f, 10*log10(event(event_sel).psd(:, channel_sel)));
plot(event(event_sel).f, event(event_sel).psd(:, channel_sel));
set(gca, 'xlim', [5 50]);
set(gca, 'ylim', [0 6]);
title(['PSD of channel ' EEG.chanlocs(channel_sel).labels ' in Attend ' stims{event_sel} ' condition'])
xlabel('Frequency (Hz)');
ylabel('Power(a.u)');



