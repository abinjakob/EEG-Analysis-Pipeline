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
LP = 45;                        % cut-off
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
EEG = pop_jointprob(EEG, 1, [1:EEG.nbchan], PRUNE, PRUNE, 0, 1, 0);
EEG = eeg_checkset(EEG);
% baseline correction
EEG = pop_rmbase(EEG, baseline);
EEG = eeg_checkset(EEG);

%% epcohing and calculating PSD for each condition

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

%% topoplots of the PSD for each events
freqs = event(1).f;
freqID(1) =  dsearchn(freqs,15);
freqID(2) =  dsearchn(freqs,20);
stims = {'Left 20Hz', 'Left 15Hz', 'Right 20Hz', 'Right 15Hz'};

% loop over events
for iEvent = 1:length(events)
    figure;
    % loop over freqID
    for idx = 1:length(freqID)
        subplot(1,2,idx);
        % extracting power to plot
        power = event(iEvent).psd(freqID(idx),:);
        % Filter out values above the upper limit
        power(power > 50) = NaN;
        % plot topoplot
        topoplot(power, EEG.chanlocs, 'electrodes', 'on', 'chaninfo', EEG.chaninfo);
        title([num2str(round(freqs(freqID(idx)))) ' Hz'])
    end 
    sgtitle(['Power for ' stims{iEvent} ' condition'])
    pause(0.1);
end 

%% plotting PSD

% params for plotting
event_sel = 1;       % event to plot
% channel_sel = [5 39 49 50];    % channel to plot
channel_sel = 23;
stims = {'Left 20Hz', 'Left 15Hz', 'Right 20Hz', 'Right 15Hz'};
% plotting PSD for selected channel
figure;
% plot(event(event_sel).f, 10*log10(event(event_sel).psd(:, channel_sel)));
plot(event(event_sel).f, event(event_sel).psd(:, channel_sel));
% plot(event(event_sel).f, mean(event(event_sel).psd(:, channel_sel), 2));
set(gca, 'xlim', [5 45]);
% set(gca, 'ylim', [0 40]);
% title(['PSD of channel ' EEG.chanlocs(channel_sel).labels ' in Attend ' stims{event_sel} ' condition'])
title(['PSD of Occipital Channels in Attend ' stims{event_sel} ' condition'])
xlabel('Frequency (Hz)');
ylabel('Power(a.u)');

%% topoplots of the attended trials 
figure;
event2plot = 4
EEG_temp = pop_selectevent(EEG, 'type', events{event2plot},'renametype', events{event2plot}, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
% plot ERP for occipital electrode
plot(EEG_temp.times,mean(mean(EEG_temp.data(23,:,:),3),1))
xlim([0 400])
figure;
pop_topoplot(EEG_temp, 1, 115, stims{event2plot}, [1 1],0, 'electrodes', 'on', 'chaninfo', EEG_temp.chaninfo);

%% weighting using GED 
% clear EEG_temp; clear narrowFilt; clear  matS; clear  matR; clear  covmatS; clear covmatR; clear evecs; clear evals; clear sidx;clear gedData;
eventSel = 3;
freqSel = 20;
% calculating for attended left 20hz 
% extracting subepoch 
EEG_temp = pop_selectevent(EEG, 'type', events{eventSel},'renametype', events{eventSel}, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');

% preparing signal matrix (matS)
% bandpass filtering data around 20Hz stim freq
narrowFilt = pop_eegfiltnew(EEG_temp, freqSel, freqSel, [], 0, [], 0);
matS = narrowFilt.data;
% reshapign data 
matS = reshape(matS, EEG_temp.nbchan, []);
% making it mean centred 
matS = bsxfun(@minus, matS, mean(matS,2));
% calculate the covariance matrix for S 
covmatS = (matS * matS')/(size(matS,2)-1);       % (dividing by a normalisation factor of n-1)

% preparing reference matrix (matR)
% broadband filtered data
matR = EEG_temp.data;
% reshapign data 
matR = reshape(matR, EEG_temp.nbchan, []);
% making it mean centred 
matR = bsxfun(@minus, matR, mean(matR,2));
% calculate the covariance matrix for S 
covmatR = (matR * matR')/(size(matR,2)-1);       % (dividing by a normalisation factor of n-1)

% generalised eigendecomposition 
[evecs, evals] = eig(covmatS, covmatR);
% sorting diagonal values of eigenvalues in ascending
[~,sidx] = sort(real(diag(evals)));
% sorting eigenvectors based on sorted eigenvalues
evecs = real(evecs(:,sidx));

% weighing the channel data based on GED
gedData = reshape( (matR'*evecs(:,end))',EEG_temp.pnts,EEG_temp.trials);

%% calculating PSD
% parameters for pwelch
window_length = size(gedData,1);
overlap = window_length / 2;

% calculating power 
% loop over trials 
for iTrial = 1:size(gedData,2)
    % computing psd usign pwelch
    [pxx, f] = pwelch(gedData(:,iTrial), hamming(window_length), overlap, 2^nextpow2(window_length*4), EEG_temp.srate);
    % computing psd usign pwelch
    pxx_all(:,iTrial) = pxx;
end

% plotting PSD
figure;
plot(f, mean(pxx_all,2));
set(gca, 'xlim', [5 45]);
stims = {'Left 20Hz', 'Left 15Hz', 'Right 20Hz', 'Right 15Hz'};
title(['PSD of GED weighted data for ' stims{eventSel} ' condition'])
xlabel('Frequency (Hz)');
ylabel('Power (norm.)');









