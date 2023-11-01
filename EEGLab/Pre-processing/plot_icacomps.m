
function plot_icacomps(DATAPATH) 
% function plot_icacomps(DATAPATH)
%
% This function is used to plot the ICA components and save it as png 
% file to the given folder.
%
% Inputs:
%   DATAPATH (char)   : function CORR_apply(DATAPATH, PATHOUT)
%
% Example function call:
% plot_icacomps(DATAPATH) 



% read all .set files in PATHIN
file_list = dir(fullfile(DATAPATH, '*.set'));

% loop over .set files 
for file_numb = 1:length(file_list)
    % extracting file names and creating subject names 
    subj{file_numb} = strrep(file_list(file_numb).name, '.set', '');
    % load rawdata to EEGLab
    EEG = pop_loadset('filename', [subj{file_numb}, '.set'], 'filepath', DATAPATH);

    % plotting ICA components for each components
    pop_topoplot(EEG, 0, [1:size(EEG.icawinv,2)], [subj{file_numb}], 0, 'electrodes', 'off');
    % creating filenames
    filename = strrep(subj{file_numb}, '_ica_comps','');
    % save as png file
    saveas(gcf, [DATAPATH, filename], 'png');
    % closing figure window
    close;
end