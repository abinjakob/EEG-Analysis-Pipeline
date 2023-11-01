
function run_ica(DATAPATH, PATHOUT, HP, HP_order, LP, LP_order, SRATE, PRUNE, PCA, PCADIMS)
% function run_ica(DATAPATH, PATHOUT, HP, HP_order, LP, LP_order, SRATE, PRUNE, PCA, PCADIMS)
%
% This function is used applying ICA in EEGLab. 
% The function loads '.set' EEG files from the given folder and prepares the
% data for ICA. First the data is filtered (low-pass & high-pass), then
% downsampled and made into epochs of 1sec. Artefactual epochs are
% eliminated and run through ICA. A PCA step is optionally performed to
% reduce ICA dimensions. After ICA the ICA weights are applied to the
% rawdata and then it is saved to a given folder. 
% 
% Inputs:
%   DATAPATH (char) : folder path for the raw data
%   PATHOUT (char)  : folder path to save the ICA filtered data
%   HP (float)      : cut-off frequency for high-pass filter 
%   HP_order (int)  : filter order for high-pass filter
%   LP (float)      : cut-off frequency for low-pass filter 
%   LP_order (int)  : filter order for low-pass filter
%   SRATE (int)     : downasampling frequency for ICA
%   PRUNE (int)     : artifact rejection threshold for artifactual epochs
%   PCA (int)       : perform PCA before ICA for dimension reduction [0 : 'No', 1: 'Yes']
%   PCADIMS (int)   : PCA dimension if PCA is opted 
%
% Example function call:
% run_ica(DATAPATH, PATHOUT, HP, HP_order, LP, LP_order, SRATE, PRUNE, PCA, PCADIMS)



% create folder if not available 
if ~exist(PATHOUT)
    mkdir(PATHOUT);
end 

% read all .set files in PATHIN
file_list = dir(fullfile(DATAPATH, '*.set'));

% loop over .set files 
for file_numb = 1:length(file_list)
    % extracting file names and creating subject names 
    subj{file_numb} = strrep(file_list(file_numb).name, '.set', '');
    % load rawdata to EEGLab
    EEG = pop_loadset('filename', [subj{file_numb}, '.set'], 'filepath', DATAPATH);
    % set name for the EEG dataset
    EEG.setname = [subj{file_numb}, '_dummy_ICA'];
    
    % prepare data for ICA
    % apply low-pass filter 
    EEG = pop_firws(EEG, 'fcutoff', LP, 'ftype', 'lowpass', 'wtype', 'hann', 'forder', LP_order);
    % downsample the data
    EEG = pop_resample(EEG, SRATE);
    % apply high-pass filter
    EEG = pop_firws(EEG, 'fcutoff', HP, 'ftype', 'highpass', 'wtype', 'hann', 'forder', HP_order);
    % create epochs for IC
    EEG = eeg_regepochs(EEG, 'recurrence', 1, 'eventtype', '999');
    EEG = eeg_checkset(EEG, 'eventconsistency');
    % remove epochs with artefacts above standard deviation PRUNE 
    EEG = pop_jointprob(EEG, 1, [1:size(EEG.data,1)], PRUNE, PRUNE, 0, 1, 0);

    % performing ICA
    % check if PCA is opted
    if PCA == 1
        % perform ICA after reducing dimension using PCA
        EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1, 'pca', PCADIMS);
    else
        % perform ICA with current dimensions
        EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1);
    end 
    
    % apply ICA weights to the data 
    % store ICA weights
    icawinv = EEG.icawinv;
    icas = EEG.icasphere;
    icaw = EEG.icaweights;
    % load the rawdata 
    EEG = pop_loadset('filename', [subj{file_numb}, '.set'], 'filepath', DATAPATH);
    
    % apply ICA weights
    EEG.icawinv = icawinv;
    EEG.icasphere = icas;
    EEG.icaweights = icaw;
    EEG = eeg_checkset(EEG);
    % set name for EEG dataset
    EEG.setname = [subj{file_numb}, '_ica'];

    % save new dataset with ICA weights
    EEG = pop_saveset(EEG, [EEG.setname, '.set'], PATHOUT);
end 