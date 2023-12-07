% % EEG DATA PRE-Processing MAIN SCRIPT 
% ------------------------------------
% The code performs the pre-processing pipeline for the artefact corrected EEG data 
% 
% Pre-processing steps:
% 1. Filtering (Low-pass and High-pass filters)
% 2. Rereference to CAR
% 3. Epoching data
% 4. Baseline Correction
% 5. Reject artefactual epochs
% 6. Calculating ERD for epoched data
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

% broad band filtering 
% high-pass filter 
HP = 1;                         % cut-off
HP_order = 826;                 % filter order    
% low-pass filter  
LP = 40;                        % cut-off
LP_order = 776;                 % filter order 

% epoching
% event markers 
events = {'left_execution','right_execution','left_imagery','right_imagery'};
epoch_start = -5;            
epoch_end = 4;               

% baseline correction
% defining baseline for baseline correcion
baseline = [epoch_start*EEG.srate 0];   

% reject artefactual epochs 
PRUNE = 4;

% parameters for ERD calculation
mu = [8 12];
beta = [13 30];
binsize = 30;
base_start = -4;
base_end = -2;

%% pre-processing 

% load ICA-clean file
%EEG = pop_loadset('L:\Cloud\NeuroCFN\RESEARCH PROJECT\EEG Analysis Scripts\temp_rawdata\MI_pilot02_ICAcleaned.set');

% low-pass filter (broad band)
EEG = pop_firws(EEG, 'fcutoff', LP, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', LP_order);
% high-pass filter (broad band)
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
EEG = pop_jointprob(EEG, 1, [1:EEG.nbchan], PRUNE, PRUNE, 0, 1, 0);
EEG = eeg_checkset(EEG);
% baseline correction
EEG = pop_rmbase(EEG, baseline);
EEG = eeg_checkset(EEG);

%% sub-epcohing and calculating ERD for each condition in Mu Band

% loop over events 
for iEvent = 1:length(events)
    % narrow-band filtering (alpha)
    EEG_temp = pop_firws(EEG, 'fcutoff', mu(2), 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', LP_order);
    EEG_temp = pop_firws(EEG_temp, 'fcutoff', mu(1), 'ftype', 'highpass', 'wtype', 'hamming', 'forder', HP_order);
    % select epoching  
    EEG_temp = pop_selectevent(EEG_temp, 'type', events{iEvent},'renametype', events{iEvent}, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
    % calculating ERD and storing
    muBand(iEvent).erd = calc_erd(EEG_temp, binsize, base_start, base_end);
    % setting event type
    muBand(iEvent).eventtype = events{iEvent};
end 

%% plotting ERD in Mu band

% creatig time points for plotting 
timevec = linspace(epoch_start,epoch_end,size(muBand(1).erd,2));
% params for plotting
event_sel = 4;               % event to plot
channel_sel = 9;             % channel to plot
% event names
stims = {'Left Execution', 'Right Execution', 'Left Imagery', 'Right Imagery'};

% plotting ERD 
figure;
plot(timevec, muBand(event_sel).erd(channel_sel,:));
title(['ERD in Mu band for channel ' EEG.chanlocs(channel_sel).labels ' in ' stims{event_sel} ' condition'])
xlabel('time (sec)');
ylabel('ERD (%)');
ylim([-100 200]);

%% sub-epcohing and calculating ERD for each condition in Beta Band

% loop over events 
for iEvent = 1:length(events)
    % narrow-band filtering (beta)
    EEG_temp = pop_firws(EEG, 'fcutoff', beta(2), 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', LP_order);
    EEG_temp = pop_firws(EEG_temp, 'fcutoff', beta(1), 'ftype', 'highpass', 'wtype', 'hamming', 'forder', HP_order);
    % select epoching  
    EEG_temp = pop_selectevent(EEG_temp, 'type', events{iEvent},'renametype', events{iEvent}, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
    % calculating ERD and storing
    betaBand(iEvent).erd = calc_erd(EEG_temp, binsize, base_start, base_end);
    % setting event type
    betaBand(iEvent).eventtype = events{iEvent};
end 

%% plotting ERD in Beta band

% creatig time points for plotting 
timevec = linspace(epoch_start,epoch_end,size(betaBand(1).erd,2));
% params for plotting
event_sel = 4;               % event to plot
channel_sel = 9;             % channel to plot
% event names
stims = {'Left Execution', 'Right Execution', 'Left Imagery', 'Right Imagery'};

% plotting ERD 
figure;
plot(timevec, betaBand(event_sel).erd(channel_sel,:));
title(['ERD in Beta band for channel ' EEG.chanlocs(channel_sel).labels ' in ' stims{event_sel} ' condition'])
xlabel('time (sec)');
ylabel('ERD (%)');
ylim([-100 200]);
