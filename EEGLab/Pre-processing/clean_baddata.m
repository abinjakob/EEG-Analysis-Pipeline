
function clean_baddata(DATAPATH, PATHOUT)
% function clean_baddata(DATAPATH, PATHOUT)
%
% This function removes the bad ICA components marked using CORRMAT
% algorithm and save the clean data
%
% Inputs:
%   DATAPATH (char)   : function CORR_apply(DATAPATH, PATHOUT)
%   PATHOUT (char)    : folder path to save the png images
%
% Example function call:
% clean_baddata(DATAPATH, PATHOUT) 



% create folder if not available 
if ~exist(PATHOUT)
    mkdir(PATHOUT);
end

% read all .set files in PATHIN
file_list = dir(fullfile(DATAPATH, '*.set'));

% loop over ICA weighted dataset
for file_numb = 1:length(file_list)
    % extracting file names and creating subject names 
    subj{file_numb} = strrep(file_list(file_numb).name, '_badcomps.set', '');
    EEG = pop_loadset('filename', [subj{file_numb}, '_badcomps.set'], 'filepath', DATAPATH);

    % remove selected componeents from dataset
    EEG = pop_subcomp(EEG, [EEG.badcomps], 0);
    EEG = eeg_checkset(EEG);

    % set name for clean EEG dataset
    EEG.setname = [subj{file_numb}, '_cleaned'];
    % save new dataset with ICA weights
    EEG = pop_saveset(EEG, 'filename', [EEG.setname, '.set'], 'filepath', PATHOUT);
end