
function convert_rawdata(DATAPATH, PATHOUT, chan_file) 
% function convert_rawdata(DATAPATH, PATHOUT, chan_file)
%
% This function is used to convert the BrainVision EEG files (.vhdr) to the
% EEGLab files (.set) for performing the data processing in EEGLab. The
% function reads all the dataset file from the given folder, converts it to
% '.set' format and saves the files to the given location
% 
% Inputs:
%   DATAPATH (char)   : folder path for the raw data
%   PATHOUT (char)    : folder path to save the converted files
%   chan_file (char)  : folder path to channle locations file
%
% Example function call:
% convert_rawdata(DATAPATH, PATHOUT, chan_file) 



% create folder if not available 
if ~exist(PATHOUT)
    mkdir(PATHOUT);
end 

% read all .vhdr files from DATAPATH
file_list = dir(fullfile(DATAPATH, '*.vhdr'));

% loop over files 
for file_numb = 1:length(file_list)
    % create a cell array of subject ids
    subj{file_numb} = ['s' sprintf('%02d', file_numb)];
    % create a cell array of names by removing .vhdr
    file_name{file_numb} = strrep(file_list(file_numb).name, '.vhdr', '');
    % load BrainVision files (.vhdr)
    EEG = pop_loadbv(DATAPATH, [file_name{file_numb}, '.vhdr']);
    % load electrode location file
    EEG = pop_chanedit(EEG, 'lookup', chan_file, 'load', {chan_file 'filetype' 'autodetect'});
    % add comments for the EEG dataset
    EEG.comments = '';
    % EEG dataset name 
    EEG.setname = [subj{file_numb}];
    % save the dataset to PATHOUT 
    EEG = pop_saveset(EEG, [EEG.setname, '.set'], PATHOUT);
end 