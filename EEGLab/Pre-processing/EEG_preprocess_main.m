% % EEG DAT PRE-Processing MAIN SCRIPT 
% ------------------------------------
% The code performs the pre-processing pipeline for the artefact corrected EEG data 
% 
% Pre-processing steps:
% 1. Filtering (Low-pass and High-pass filters)
% 2. Rereference to CAR
% 3. Epoching data
% 4. Baseline Correction
% 5. Reject artefactual epochs
% 
% Author: Abin Jacob
% Date  : 04/11/2023

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
addpath('/Users/abinjacob/Documents/02. NeuroCFN/eeglab2023.1');
% open EEGLab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% set current directory 
cd(MAINPATH);

%% SET PARAMETERS FOR PRE-PROCESSING

% filtering 
% high-pass filter 
HP = 0.1;                       % cut-off
HP_order = 500;                 % filter order    
% low-pass filter  
LP = 40;                        % cut-off
LP_order = 100;                 % filter order 

% epoching 
event_name = 'S 10';            % event name for the stim onset   
epoch_start = -0.2;             % epoch start  
epoch_end = 0.8;                % epoch end 

% baseline correction
fs = 1000;                      % sampling rate of the EEG recording 
baseline = [epoch_start*fs 0];    % defining baseline for the correction

% reject artefactual epochs 
PRUNE = 4;

%% pre-processing for all subjects

% create folder if not available 
if ~exist(PATHOUT)
    mkdir(PATHOUT);
end

% read all .set files in PATHIN
file_list = dir(fullfile(DATAPATH, '*.set'));

% loop over ICA weighted dataset
for file_numb = 1:length(file_list)
    % extracting file names and creating subject names 
    subj{file_numb} = strrep(file_list(file_numb).name, '_ica_cleaned.set', '');
    EEG = pop_loadset('filename', [subj{file_numb}, '_ica_cleaned.set'], 'filepath', DATAPATH);

    % filtering 
    % low-pass filter 
    EEG = pop_firws(EEG, 'fcutoff', LP, 'ftype', 'lowpass', 'wtype', 'hann', 'forder', LP_order);
    % high-pass filter 
    EEG = pop_firws(EEG, 'fcutoff', HP, 'ftype', 'highpass', 'wtype', 'hann', 'forder', HP_order);

    % re-referencing to CAR
    EEG = pop_reref(EEG, [], 'refstate',0);

    % removing unnecessary event marker
    event_pos = 1;      % position counter for the events other than stim onset
    event_idx = [];     % array to store the index of the event other than stim onset

    for idx = 1: length(EEG.event)
        if ~ strcmp(EEG.event(idx).type, event_name)
            event_idx(event_pos) = idx;
            event_pos = event_pos + 1;
        end
    end 
    % remove events which are not stim onset from the data
    EEG = pop_editeventvals(EEG, 'delete', event_idx);
    EEG = eeg_checkset(EEG);

    % epoching 
    EEG = pop_epoch(EEG, {event_name}, [epoch_start epoch_end], 'newname', EEG.setname, 'epochinfo', 'yes');
    % baseline correction
    EEG = pop_rmbase(EEG, baseline);

    % reject artefactual epochs 
    % joint probability-based artifact rejection (joint prob. > PRUNE (SD))
    EEG = pop_jointprob(EEG, 1, [1:EEG.nbchan], PRUNE, PRUNE, 0, 0, 0);
    % update EEG data srtucture based on artefact rejection
    EEG = eeg_rejsuperpose(EEG, 1, 1, 1, 1, 1, 1, 1, 1);
    % store rejected epochs 
    rej_epoch (file_numb) = {find(EEG.reject.rejglobal == 1)};
    % reject bad epochs
    EEG = pop_rejepoch(EEG, EEG.reject.rejglobal, 0);

    % get indices of ICA components that were rejected 
    ica_rejcomps(file_numb) = {EEG.badcomps};

    % save pruned and epoched dataset 
    EEG = eeg_checkset(EEG);
    % set name for the EEG dataset
    EEG.setname = [subj{file_numb}, '_epoched'];
    EEG = pop_saveset(EEG, [EEG.setname, '.set'], PATHOUT);
end


