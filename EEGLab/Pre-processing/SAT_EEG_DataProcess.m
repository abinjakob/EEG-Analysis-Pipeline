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
% 6. Sub epoch the data to attend left and right conditions
% 7. Calculate ERP
% 
% the pre-processed data is then stored to the given folder
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
addpath('L:\Cloud\SW\eeglab2023.1');
% open EEGLab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% set current directory 
cd(MAINPATH);

%% SET PARAMETERS FOR PRE-PROCESSING

% filtering 
% high-pass filter 
HP = 0.1;                       % cut-off
HP_order = 826;                 % filter order    
% low-pass filter  
LP = 40;                        % cut-off
LP_order = 776;                 % filter order 

% epoching 
% extracting correct events
[trialCount, corrEvent_left, corrEvent_right] = SATcorrect_trials(EEG);
event_name = [corrEvent_left, corrEvent_right];
epoch_start = -0.25;             % epoch start  
epoch_end = 3.25;                  % epoch end 

% baseline correction
fs = EEG.srate;                 % sampling rate of the EEG recording 
baseline = [epoch_start*fs 0];  % defining baseline for the correction

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
    EEG = pop_epoch(EEG, event_name, [epoch_start epoch_end], 'newname', 'SET_pilot_Epoched', 'epochinfo', 'yes');

    % reject artefactual epochs 
    % joint probability-based artifact rejection (joint prob. > PRUNE (SD))
    EEG = pop_jointprob(EEG, 1, [1:EEG.nbchan], PRUNE, PRUNE, 0, 1, 0);
    % save pruned and epoched dataset 
    EEG = eeg_checkset(EEG);

    % baseline correction
    EEG = pop_rmbase(EEG, baseline);
    % save pruned and epoched dataset 
    EEG = eeg_checkset(EEG);
    % set name for the EEG dataset
    EEG.setname = [subj{file_numb}, '_epoched'];
    EEG = pop_saveset(EEG, [EEG.setname, '.set'], PATHOUT);
end

%% calculating ERP for left and right attended trials 

% parameters for calculating ERP
% baseline period 
baseline_period = size(baseline(1):baseline(2)-1,2);
% stimulus period 
stim_duration = 3;
% number of left and right tones
left_toneCount = 4;
right_toneCount = 5;
% tone duration
left_toneDuration = (stim_duration*EEG.srate)/left_toneCount;
right_toneDuration = (stim_duration*EEG.srate)/right_toneCount;
% 
toneStart = baseline_period + 1;

% selecting left attended epochs 
EEG_left = pop_selectevent(EEG, 'type',corrEvent_left,'renametype','left','deleteevents','off','deleteepochs','on','invertepochs','off');
% selecting right attended epochs 
EEG_right = pop_selectevent(EEG, 'type',corrEvent_right,'renametype','right','deleteevents','off','deleteepochs','on','invertepochs','off');

%% for attend left condition

% calculate ERP for left tones 
% loop over channels 
for iChannel = 1:size(EEG_left.data,1)
    for iSegment = 1:left_toneCount-1
        start = (toneStart + (left_toneDuration * iSegment));
        stop = (start + left_toneDuration) - 1;
        attendLeft_lefttone(iChannel,:,iSegment) = squeeze(mean(EEG_left.data(iChannel,start:stop,:), 3)); 
    end 
end 

% calculate ERP for right tones
% loop over channels 
for iChannel = 1:size(EEG_left.data,1)
    for iSegment = 1:right_toneCount-1
        start = (toneStart + (right_toneDuration * iSegment));
        stop = (start + right_toneDuration) - 1;
        attendLeft_righttone(iChannel,:,iSegment) = squeeze(mean(EEG_left.data(iChannel,start:stop,:), 3)); 
    end 
end 

% plotting the ERP
figure;
% plotting the left tones
plot(EEG_left.times(toneStart:(toneStart+left_toneDuration)-1),mean(mean(attendLeft_lefttone([2 5 6 7 8 9 ],:,:),3),1))
hold on
% plotting the right tones
plot(EEG_left.times(toneStart:(toneStart+right_toneDuration)-1),mean(mean(attendLeft_righttone([2 5 6 7 8 9 ],:,:),3),1))
hold off
% setting lims 
xlim([0 600]);
ylim([-1.6 0.2]);
xlabel('time (ms)');
ylabel('amplitude (mV)');
title('ERP for central electrodes for Attend Left Condition');
legend('Left Tones', 'Right Tones');

% 
% figure;
% % plotting the left tones
% plot(EEG_left.times(toneStart:(toneStart+left_toneDuration)-1),mean(attendLeft_lefttone([2 5 6 7 8 9 ],:,3),1))
% hold on
% % plotting the right tones
% plot(EEG_left.times(toneStart:(toneStart+right_toneDuration)-1),mean(attendLeft_righttone([2 5 6 7 8 9 ],:,3),1))
% hold off
% % setting lims 
% xlim([0 600]);
% % ylim([-1.4 0.2]);

%% for attend right condition

% calculate ERP for left tones 
% loop over channels 
for iChannel = 1:size(EEG_right.data,1)
    for iSegment = 1:left_toneCount-1
        start = (toneStart + (left_toneDuration * iSegment));
        stop = (start + left_toneDuration) - 1;
        attendRight_lefttone(iChannel,:,iSegment) = squeeze(mean(EEG_right.data(iChannel,start:stop,:), 3)); 
    end 
end 

% calculate ERP for right tones
% loop over channels 
for iChannel = 1:size(EEG_right.data,1)
    for iSegment = 1:right_toneCount-1
        start = (toneStart + (right_toneDuration * iSegment));
        stop = (start + right_toneDuration) - 1;
        attendRight_righttone(iChannel,:,iSegment) = squeeze(mean(EEG_right.data(iChannel,start:stop,:), 3)); 
    end 
end 

% plotting the ERP
figure;
% plotting the left tones
plot(EEG_right.times(toneStart:(toneStart+left_toneDuration)-1),mean(mean(attendRight_lefttone([2 5 6 7 8 9 ],:,:),3),1))
hold on
% plotting the right tones
plot(EEG_right.times(toneStart:(toneStart+right_toneDuration)-1),mean(mean(attendRight_righttone([2 5 6 7 8 9 ],:,:),3),1))
hold off
% setting xlims 
xlim([0 600]);
ylim([-1.6 0.2]);
xlabel('time (ms)');
ylabel('amplitude (mV)');
title('ERP for central electrodes for Attend Right Condition');
legend('Left Tones', 'Right Tones');

% figure; plot(EEG_right.times,mean(mean(EEG_right.data([2 5 6 7 8 9 ],:,:),3),1))
% xlim([0 3000]); ylim([-3 1.5]);

%% plot central electrodes 
figure; plot(EEG_left.times,mean(mean(EEG_left.data([2 5 6 7 8 9 ],:,:),3),1))
hold on
plot(EEG_right.times,mean(mean(EEG_right.data([2 5 6 7 8 9 ],:,:),3),1))
plot(0:3000/4:3000,[0 0 0 0 0],'db')
plot(0:3000/5:3000,[0 0 0 0 0 0],'dr')
xlim([-250 3250]); ylim([-3 1.5]);
xlabel('time (ms)');
ylabel('amplitude (mV)');
title('ERP for central electrodes for Attend Left and Attend Right Condition');
legend('Attend Left', 'Attend Right', 'Left Tone Onsets', 'Right Tone Onsets');
hold off

%% plot all electrodes for attend left

figure; plot(EEG_left.times,mean(EEG_left.data(:,:,:),3))
hold on
plot(0:3000/4:3000,[-6 -6 -6 -6 -6],'db')
plot(0:3000/5:3000,[-6 -6 -6 -6 -6 -6],'dr')
xlim([-250 3250]); ylim([-6 8]);
xlabel('time (ms)');
ylabel('amplitude (mV)');
title('ERP of all electrodes for Attend Left Condition');
hold off

%% plot all electrodes for attend right

figure; plot(EEG_left.times,mean(EEG_right.data(:,:,:),3))
hold on
plot(0:3000/4:3000,[-6 -6 -6 -6 -6],'db')
plot(0:3000/5:3000,[-6 -6 -6 -6 -6 -6],'dr')
xlim([-250 3250]); ylim([-6 8]);
xlabel('time (ms)');
ylabel('amplitude (mV)');
title('ERP of all electrodes for Attend Right Condition');
hold off
