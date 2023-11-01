
function CORR_apply(DATAPATH, PATHOUT)
% function CORR_apply(DATAPATH, PATHOUT)
%
% This function applies the artefactual components identified using CORRMAP 
% algorithm to the iCA weighted data marks them as bad components in the 
% EEGLab structure. The function also plots the ICA components along with
% the bad components listed and save the figure as png.
%
% Inputs:
%   DATAPATH (char)   : function CORR_apply(DATAPATH, PATHOUT)
%   PATHOUT (char)    : folder path to save the png images
%
% Example function call:
% CORR_apply(DATAPATH, PATHOUT) 



% create folder if not available 
if ~exist(PATHOUT)
    mkdir(PATHOUT);
end

% load components generated using CORRMAP algorithm
load(fullfile(DATAPATH, 'components.mat'));

% read all .set files in PATHIN
file_list = dir(fullfile(DATAPATH, '*.set'));

% loop over ICA weighted dataset
for file_numb = 1:length(file_list)
    % extracting file names and creating subject names 
    subj{file_numb} = strrep(file_list(file_numb).name, '.set', '');
    EEG = pop_loadset('filename', [subj{file_numb}, '.set'], 'filepath', DATAPATH);

    % find selected artefact components for the dataset
    % eye blink
    eye = eyeblink(find(eyeblink(:,1) == file_numb),2)';
    % heart beat
    heart = heartbeat(find(heartbeat(:,1) == file_numb),2)';
    % lateral eye movement
    eyemov = eyemovement(find(eyemovement(:,1) == file_numb),2)';

    % store all components in EEG structure as bad components
    EEG.badcomps = [eye heart eyemov];
    
    % set name for EEG dataset
    EEG.setname = [subj{file_numb}, '_badcomps'];
    % save new dataset with ICA weights
    EEG = pop_saveset(EEG, 'filename', [EEG.setname, '.set'], 'filepath', PATHOUT);

    % plot the components with bad componenets labelled
    pop_topoplot(EEG, 0, [1:size(EEG.icawinv, 2)], [subj{file_numb}, ' removed components: ', num2str(EEG.badcomps)], ...
        [6,10], 0, 'electrodes', 'off');
    % filename for saving figure 
    filename = strrep(subj{file_numb}, '_ica', '_badcomps');
    % save figure as png
    saveas(gcf, [PATHOUT, filename], 'png');
    close;
end