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
epoch_end = 6;
% epoch_end = 4;

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
event_sel = 3;               % event to plot
channel_sel = 9;             % channel to plot
% event names
stims = {'Left Execution', 'Right Execution', 'Left Imagery', 'Right Imagery'};

% plotting ERD 
figure;
plot(timevec, muBand(event_sel).erd(channel_sel,:));
title(['ERD in Mu band for channel ' EEG.chanlocs(channel_sel).labels ' in ' stims{event_sel} ' condition'])
xlabel('time (sec)');
ylabel('ERD (%)');
xlim([-5 6])
ylim([-100 200]);
rectangle('Position', [0, -100, 4, 300], 'EdgeColor', 'none', 'FaceColor', [1 0 0 0.03]);

%% plotting ERD for left and right electrodes (for high density cap)

% creatig time points for plotting 
timevec = linspace(epoch_start,epoch_end,size(muBand(1).erd,2));
% event names
stims = {'Left Execution', 'Right Execution', 'Left Imagery', 'Right Imagery'};
% channels to plot (left and right electrodes near M1 area)
chan2plot = [35 36 7 41 6 56; 32 33 38 4 3 48];

% loop over events
for iEvent = 1:length(events)
    % opening figure window for each event
    figure('Position', [10 10 900 600]);
    % loop over left and right 
    for iPos = 1:size(chan2plot,1)
        % choose subplot position
        subplot(1,2,iPos);
        % loop over channels
        for iChan = 1:size(chan2plot,2)
            % choose channel to plot
            chanVal = chan2plot(iPos,iChan);
            % plottting ERD of different channels in same plot
            plot(timevec, muBand(iEvent).erd(chanVal,:));
            hold on
        end 
        % set title for subplot
        if iPos == 1
            title('Left Electrodes near M1')
        else
            title('Right Electrodes near M1')
        end 
        % set plot parameters
        xlabel('time (sec)');
        ylabel('ERD (%)');
        xlim([-5 6])
        ylim([-100 200]);
        rectangle('Position', [0, -100, 4, 300], 'EdgeColor', 'none', 'FaceColor', [1 0 0 0.03]);
        % stop plotting in the same plot and moving to next subplot
        hold off
    end 
    % set title for the figure
    sgtitle(['Mu band ERD in ' stims{iEvent} ' condition'])
end 

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
event_sel = 1;               % event to plot
channel_sel = 35;             % channel to plot
% event names
stims = {'Left Execution', 'Right Execution', 'Left Imagery', 'Right Imagery'};

% plotting ERD 
figure;
plot(timevec, betaBand(event_sel).erd(channel_sel,:));
title(['ERD in Beta band for channel ' EEG.chanlocs(channel_sel).labels ' in ' stims{event_sel} ' condition'])
xlabel('time (sec)');
ylabel('ERD (%)');
ylim([-100 200]);

%% topoplots 

bin2plot = 120;
event2plot = 4;
% plotting topoplot for mu band
figure; topoplot(muBand(event2plot).erd(:,bin2plot), EEG.chanlocs, 'electrodes', 'on', 'chaninfo', EEG.chaninfo);


%% ------------------------------------------------
%% weigting using GED 

% event to calculate
event_sel = 3;
% pre-stim period
pre_dur = [-5 -1];              % in sec
pre_dur = pre_dur * 1000;       % converted to milli sec 
% stim period 
stim_dur = [0 4];               % in sec
stim_dur = stim_dur * 1000;     % converted to milli sec 
% index of the stim and pre-stim periods
pre_Tidx = dsearchn(EEG.times',[pre_dur(1) pre_dur(2)]');  
stim_Tidx = dsearchn(EEG.times',[stim_dur(1) stim_dur(2)]');  

% narrowband filtering (mu)
EEG_temp = pop_firws(EEG, 'fcutoff', mu(2), 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', LP_order);
EEG_temp = pop_firws(EEG_temp, 'fcutoff', mu(1), 'ftype', 'highpass', 'wtype', 'hamming', 'forder', HP_order);
% select event epochs 
EEG_temp = pop_selectevent(EEG_temp, 'type', events{event_sel},'renametype', events{event_sel}, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');

% preparing signal matrix (matS)
matS = reshape(EEG_temp.data(:,stim_Tidx(1):stim_Tidx(2),:),EEG_temp.nbchan,[]);
% making it mean centred 
matS = bsxfun(@minus, matS, mean(matS,2));
% calculate the covariance matrix for S 
covmatS = (matS * matS')/(size(matS,2)-1);       % (dividing by a normalisation factor of n-1)

% preparing reference matrix (matS)
matR = reshape(EEG_temp.data(:,pre_Tidx(1):pre_Tidx(2),:),EEG_temp.nbchan,[]);
% making it mean centred 
matR = bsxfun(@minus, matR, mean(matR,2));
% calculate the covariance matrix for S 
covmatR = (matR * matR')/(size(matR,2)-1);       % (dividing by a normalisation factor of n-1)

% generalised eigendecomposition 
[evecs, evals] = eig(covmatS, covmatR);
% sorting diagonal values of eigenvalues in ascending
[~,eigidx] = sort(real(diag(evals)));
% sorting eigenvectors based on sorted eigenvalues
evecs = real(evecs(:,eigidx));

% weighing the channel data based on GED
gedData = reshape(EEG_temp.data, EEG_temp.nbchan, [])' * evecs(:,end);
gedData = reshape(gedData, EEG_temp.pnts, EEG_temp.trials);

%% Calculate ERD for GED weighted channels

% window size
window = EEG_temp.pnts/binsize; 
% calculating time points 
tbins = [];
% loop over window
for t = 1:window
    start = 1+(t-1)*binsize;
    tbins{t} = (start:(start+binsize)-1);
end 

% squaring the amplitude to obtain power (step 02)
% loop over bins 
for iBin = 1:length(tbins)    
    % loop over trials
    for iTrial = 1:size(gedData,2)
        chanPower_vals(iBin,iTrial,:) = gedData(tbins{iBin},iTrial).^2;
    end
end

% averaging across trials (step 03)
chanPower_TrialsAvg = squeeze(mean(chanPower_vals,2));
% averaging across time samples (step 04) (period after event, A)
chanPower = mean(chanPower_TrialsAvg,2);

% baseline duration
baseline_period = base_end - base_start;

% calculate ERD
% calculating avg power in baseline period (reference period, R)
baseline_avg = mean(chanPower(1:(floor((baseline_period*EEG.srate)/binsize))));
% calculating ERD% = ((A-R)/R)*100
erd = ((chanPower-baseline_avg)/baseline_avg)*100;

%% plotting ERD after GED

% creatig time points for plotting 
timevec = linspace(epoch_start,epoch_end,size(erd,1));
% event names
stims = {'Left Execution', 'Right Execution', 'Left Imagery', 'Right Imagery'};

% plotting ERD 
figure;
plot(timevec, erd);
title(['ERD in Mu band after GED in ' stims{event_sel} ' condition'])
xlabel('time (sec)');
ylabel('ERD (%)');
xlim([-5 4])
ylim([-100 200]);
% rectangle('Position', [0, -100, 4, 300], 'EdgeColor', 'none', 'FaceColor', [1 0 0 0.03]);

%% plotting topoplot
maps = inv(evecs');
figure;
topoplot(maps(:,eigidx(end)),EEG_temp.chanlocs,'numcontour',0);


