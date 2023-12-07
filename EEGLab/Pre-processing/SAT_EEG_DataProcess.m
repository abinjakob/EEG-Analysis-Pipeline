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
event_name = {'left_alt_116'  'left_alt_118'  'left_alt_12'  'left_alt_127'  'left_alt_170'  'left_alt_171'  'left_alt_183'  'left_alt_189'  'left_alt_190'  'left_alt_23'  'left_alt_280'  'left_alt_281'  'left_alt_289'  'left_alt_291'  'left_alt_294'  'left_alt_298'  'left_alt_332'  'left_alt_343'  'left_alt_385'  'left_alt_443'  'left_alt_444'  'left_alt_79'  'left_alt_8'  'left_alt_80'  'left_alt_83'  'left_asc_123'  'left_asc_143'  'left_asc_16'  'left_asc_166'  'left_asc_18'  'left_asc_222'  'left_asc_252'  'left_asc_271'  'left_asc_302'  'left_asc_303'  'left_asc_304'  'left_asc_340'  'left_asc_357'  'left_asc_358'  'left_asc_36'  'left_asc_393'  'left_asc_433'  'left_asc_434'  'left_asc_448'  'left_asc_463'  'left_asc_465'  'left_asc_68'  'left_dec_102'  'left_dec_107'  'left_dec_157'  'left_dec_159'  'left_dec_162'  'left_dec_199'  'left_dec_209'  'left_dec_258'  'left_dec_259'  'left_dec_260'  'left_dec_265'  'left_dec_269'  'left_dec_315'  'left_dec_320'  'left_dec_321'  'left_dec_324'  'left_dec_362'  'left_dec_366'  'left_dec_370'  'left_dec_371'  'left_dec_378'  'left_dec_419'  'left_dec_46'  'left_dec_471'  'left_dec_472'  'left_dec_474'  'left_dec_51'  'left_dec_93'  'right_alt_167'  'right_alt_168'  'right_alt_172'  'right_alt_178'  'right_alt_193'  'right_alt_202'  'right_alt_204'  'right_alt_206'  'right_alt_207'  'right_alt_208'  'right_alt_215'  'right_alt_219'  'right_alt_234'  'right_alt_237'  'right_alt_240'  'right_alt_244'  'right_alt_247'  'right_alt_266'  'right_alt_268'  'right_alt_57'  'right_alt_60'  'right_alt_62'  'right_alt_65'  'right_alt_72'  'right_alt_73'  'right_alt_74'  'right_alt_90'  'right_alt_98'  'right_alt_99'  'right_asc_109'  'right_asc_121'  'right_asc_125'  'right_asc_14'  'right_asc_140'  'right_asc_151'  'right_asc_153'  'right_asc_156'  'right_asc_17'  'right_asc_275'  'right_asc_283'  'right_asc_305'  'right_asc_309'  'right_asc_311'  'right_asc_313'  'right_asc_318'  'right_asc_33'  'right_asc_4'  'right_asc_47'  'right_asc_49'  'right_asc_5'  'right_dec_326'  'right_dec_327'  'right_dec_333'  'right_dec_336'  'right_dec_337'  'right_dec_354'  'right_dec_355'  'right_dec_359'  'right_dec_367'  'right_dec_372'  'right_dec_374'  'right_dec_390'  'right_dec_395'  'right_dec_398'  'right_dec_399'  'right_dec_406'  'right_dec_415'  'right_dec_426'  'right_dec_430'  'right_dec_442'  'right_dec_449'  'right_dec_455'  'right_dec_473'  'right_dec_476'  'right_dec_481'};
%event_name = {'left_alt_116'  'left_alt_127'  'left_alt_133'  'left_alt_135'  'left_alt_19'  'left_alt_190'  'left_alt_22'  'left_alt_225'  'left_alt_239'  'left_alt_241'  'left_alt_245'  'left_alt_289'  'left_alt_291'  'left_alt_296'  'left_alt_30'  'left_alt_345'  'left_alt_348'  'left_alt_351'  'left_alt_399'  'left_alt_442'  'left_alt_444'  'left_alt_451'  'left_alt_61'  'left_alt_7'  'left_asc_1'  'left_asc_114'  'left_asc_123'  'left_asc_125'  'left_asc_175'  'left_asc_194'  'left_asc_218'  'left_asc_222'  'left_asc_232'  'left_asc_250'  'left_asc_251'  'left_asc_273'  'left_asc_275'  'left_asc_305'  'left_asc_31'  'left_asc_328'  'left_asc_33'  'left_asc_338'  'left_asc_358'  'left_asc_382'  'left_asc_384'  'left_asc_392'  'left_asc_411'  'left_asc_448'  'left_asc_449'  'left_asc_6'  'left_dec_101'  'left_dec_149'  'left_dec_153'  'left_dec_157'  'left_dec_200'  'left_dec_203'  'left_dec_205'  'left_dec_209'  'left_dec_259'  'left_dec_261'  'left_dec_264'  'left_dec_265'  'left_dec_315'  'left_dec_317'  'left_dec_361'  'left_dec_365'  'left_dec_415'  'left_dec_416'  'left_dec_422'  'left_dec_423'  'left_dec_424'  'left_dec_425'  'left_dec_47'  'left_dec_480'  'left_dec_54'  'right_alt_100'  'right_alt_104'  'right_alt_168'  'right_alt_170'  'right_alt_176'  'right_alt_183'  'right_alt_184'  'right_alt_185'  'right_alt_186'  'right_alt_189'  'right_alt_191'  'right_alt_196'  'right_alt_206'  'right_alt_210'  'right_alt_212'  'right_alt_215'  'right_alt_216'  'right_alt_223'  'right_alt_224'  'right_alt_226'  'right_alt_227'  'right_alt_243'  'right_alt_244'  'right_alt_254'  'right_alt_270'  'right_alt_64'  'right_alt_66'  'right_alt_74'  'right_alt_76'  'right_alt_78'  'right_alt_82'  'right_alt_84'  'right_alt_92'  'right_alt_95'  'right_alt_99'  'right_asc_117'  'right_asc_121'  'right_asc_124'  'right_asc_126'  'right_asc_128'  'right_asc_13'  'right_asc_131'  'right_asc_134'  'right_asc_138'  'right_asc_139'  'right_asc_145'  'right_asc_159'  'right_asc_161'  'right_asc_17'  'right_asc_272'  'right_asc_284'  'right_asc_295'  'right_asc_304'  'right_asc_32'  'right_asc_36'  'right_dec_337'  'right_dec_339'  'right_dec_353'  'right_dec_355'  'right_dec_362'  'right_dec_369'  'right_dec_385'  'right_dec_400'  'right_dec_408'  'right_dec_413'  'right_dec_434'  'right_dec_435'  'right_dec_440'  'right_dec_445'  'right_dec_454'  'right_dec_462'  'right_dec_468'  'right_dec_474'  'right_dec_478'  'right_dec_479'};
epoch_start = -0.2;             % epoch start  
epoch_end = 3;                  % epoch end 

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
stim_duration = epoch_end;
% number of left and right tones
left_toneCount = 4;
right_toneCount = 5;
% tone duration
left_toneDuration = (stim_duration*EEG.srate)/left_toneCount;
right_toneDuration = (stim_duration*EEG.srate)/right_toneCount;
% 
toneStart = baseline_period + 1;

% selecting left attended epochs 
EEG_left = pop_selectevent(EEG, 'type',{event_name{1:75}},'renametype','left','deleteevents','off','deleteepochs','on','invertepochs','off');
% selecting right attended epochs 
EEG_right = pop_selectevent(EEG, 'type',{event_name{76:end}},'renametype','right','deleteevents','off','deleteepochs','on','invertepochs','off');

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
ylim([-1.4 0.2]);
xlabel('time (ms)');
ylabel('amplitude (mV)');
title('ERP for central electrodes for Attend Left Condition');
legend('Left Tones', 'Right Tones');

% % TRASH
% figure; plot(EEG_left.times,mean(mean(EEG_left.data([2 5 6 7 8 9 ],:,:),3),1))
% xlim([0 600]);
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
plot(EEG_right.times(toneStart:(toneStart+left_toneDuration)-1),mean(mean(attendRight_lefttone(5:9,:,:),3),1))
hold on
% plotting the right tones
plot(EEG_right.times(toneStart:(toneStart+right_toneDuration)-1),mean(mean(attendRight_righttone(5:9,:,:),3),1))
hold off
% setting xlims 
xlim([0 600]);
ylim([-1.4 0.2]);
xlabel('time (ms)');
ylabel('amplitude (mV)');
title('ERP for central electrodes for Attend Right Condition');
legend('Left Tones', 'Right Tones');

% figure; plot(EEG_right.times,mean(mean(EEG_right.data([2 5 6 7 8 9 ],:,:),3),1))
% xlim([0 600]);


