% EEG ICA ARTEFACT ATTENUATION MAIN SCRIPT 
% ----------------------------------------
% The code performs the artefact attenuation using ICA. 
% 
% It uses various custom functions to first convert the BrainVision data file 
% to a EEGLab readable format and prepares the data to perform ICA and use 
% CORRMAP algorithm to label the ICA components as artefacts.
% 
% Following functions are used in the script:
% (Use help to know how to use the functions)
% - convert_rawdata : convert BrainVision (.vhdr) to EEGLab files structure (.set)
% - run_ica         : run ICA on a dummy data and then applies the weights to the original dataset
% - plot_icacomps   : plots topographic map of the ICA components
% - CORR_run        : runs CORRMAP algorithm to find correlation between ICA components and the selected template
% - CORR_apply      : applies the CORRMAP values on the dataset and plots the bad componenets
% - clean_baddata   : removes the badcomponents and saves the cleaned data
% 
% Author: Abin Jacob
% Date  : 29/10/2023

%% start fresh 
clear; clc; close all;

%% define project folders

% path to the main folder 
MAINPATH = '/Users/abinjacob/Documents/02. NeuroCFN/Research Module/Paradigm/Analysis Codes/EEG Data Analysis';
% path for the rawdata
DATAPATH = fullfile(MAINPATH,'rawdata_vhdr', filesep);
% path to save output
PATHOUT = fullfile(MAINPATH,'rawdata', filesep);
% file for channel location
chan_file = fullfile(MAINPATH,'config', 'elec_64ch.elp');

% add EEGLab to matlab path
addpath('L:\Cloud\SW\eeglab2023.1');
% open EEGLab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%% convert BrainVision files to EEGLab files
cd(MAINPATH);
convert_rawdata(DATAPATH, PATHOUT, chan_file);

%% convert XDF files to EEGLab files 
cd(MAINPATH);
chan_file = fullfile(MAINPATH,'config', 'mobile24.elp');
convert_rawxdf(DATAPATH, PATHOUT, chan_file) 

% remove any faulty channel (optional)
EEG = pop_select( EEG, 'rmchannel',{'Fp2'});

%% ICA BASED ARTEFACT ATTENUATION 
%% run ICA and apply weights to rawdata

% set input and output paths
DATAPATH = fullfile(MAINPATH,'rawdata', filesep);;
PATHOUT = fullfile(MAINPATH, 'ica_corrected_data','ICA_weighted',filesep);

% set parameters for ICA
% high-pass filter 
HP = 1;                   % cut-off
HP_order = 776;           % filter order    
% low-pass filter  
LP = 40;                  % cut-off
LP_order = 776;           % filter order 
% downsampling freq. for ICA
SRATE = 250;
% artifact rejection threshold on SD
PRUNE = 3;
% perform PCA before ICA for dimension reduction [0 : 'No', 1: 'Yes']
PCA = 0;
% PCA dimension 
PCADIMS = 50;

% running ICA and applying weights to the rawdata
cd(MAINPATH);
run_ica(DATAPATH, PATHOUT, HP, HP_order, LP, LP_order, SRATE, PRUNE, PCA, PCADIMS)

%% setting paths
DATAPATH = fullfile(MAINPATH, 'ica_corrected_data','ICA_weighted',filesep);
PATHOUT = fullfile(MAINPATH, 'ica_corrected_data','ICA_badcomps',filesep);

%% plotting ICA components
plot_icacomps(DATAPATH)

%% creating STUDY for CORRMAP

% study name
studyname = 'ICA_Comps';
% initialise study index 
idx = 1;
% reset EEGLab variables
STUDY = [];
CURRENTSTUDY = 0;
ALLEEG = [];
EEG = [];
CURRENTSET = [];

% allowing to process more than one dataset at a time
pop_editoptions('option_storedisk',1);

% read all .set files in PATHIN
file_list = dir(fullfile(DATAPATH, '*.set'));

% loop over .set files 
for file_numb = 1:length(file_list)
    % extracting file names and creating subject names 
    subj{file_numb} = strrep(file_list(file_numb).name, '.set', '');
    dataset = [DATAPATH, subj{file_numb}, '.set'];
    [STUDY ALLEEG] = std_editset(STUDY, ALLEEG, 'name', studyname, 'commands', {{'index' idx 'load' dataset 'subject' subj{file_numb}}}, ...
        'updatedat', 'off', 'savedat', 'off', 'filename', [DATAPATH, studyname]);
    % increment index
    idx = idx+1;
end 

CURRENTSTUDY = 1;
EEG = ALLEEG;
CURRENTSET = [1:length(EEG)];
STUDY.design = [];
[STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG)
eeglab redraw

%% NOTE:
% choose the template topography from the ICA component plot for CORRMAP
% each time ICA is run the components order might chnage. So make sure to
% give the right components for the CORRMAP algorithm below
% CORRMAP clusters the ICA components for each dataset based on a template topography that was 
% manually selected. Choose a template topography manually from the dataset topographic ICA component 
% plots generated using 'plot_icacomps' fn and mark eyeblink, heart beat and lateral eye movement 
% components in it. CORRMAP will then compares this template with other dataset and marks the 
% topography with similar artefact components.
%% find ICA components using CORRMAP

% setting parameters for CORRMAP template
% eye blink component
blink_set = 3;              % template set choosen
blink_comp = 1;             % component for eye blink in the template set
% heart beat component
heart_set = 3;              % template set choosen
heart_comp = 9;            % component for heart beat in the template set
% eye movement component
eyemov_set = 3;             % template set choosen
eyemov_comp = 22;           % component for lateral eye movement in the template set

%% running CORRMAP and applying to data

% CORR_run also saves the components identified to a MAT-File
CORR_run(DATAPATH, studyname, blink_set, blink_comp, heart_set, heart_comp,eyemov_set, eyemov_comp)
% apply artefact componenets to ica weighted data and plot 
CORR_apply(DATAPATH, PATHOUT)

%% reject badcomponets from data

% setting paths
DATAPATH = fullfile(MAINPATH, 'ica_corrected_data','ICA_badcomps',filesep);;
PATHOUT = fullfile(MAINPATH, 'ica_corrected_data','ICA_cleaned',filesep);
% remove artefacts
clean_baddata(DATAPATH, PATHOUT);



