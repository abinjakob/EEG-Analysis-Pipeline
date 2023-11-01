
function CORR_run(DATAPATH, studyname, blink_set, blink_comp, heart_set, heart_comp,eyemov_set, eyemov_comp)
% function CORR_run(DATAPATH, studyname, blink_set, blink_comp, heart_set, heart_comp,eyemov_set, eyemov_comp)
%
% This function is used run CORRMAP alogorithm on the dataset to identify
% artefact components from IC components.
%
% CORRMAP clusters the ICA components for each dataset based on a template topography that was 
% manually selected. Choose a template topography manually from the dataset topographic ICA component 
% plots generated using 'plot_icacomps' fn and mark eyeblink, heart beat and lateral eye movement 
% components in it. CORRMAP will then compares this template with other dataset and marks the 
% topography with similar artefact components.
% 
% Inputs:
%   DATAPATH (char)   : function CORR_apply(DATAPATH, PATHOUT)
%   studyname (char)  : name of the STUDY
%   blink_set (int)   : template file choosen for eye blink
%   blink_comp (int)  : component for eye blink in the template file
%   heart_set (int)   : template file choosen for heart rate
%   heart_comp (int)  : component for lateral heart rate in the template file
%   eyemov_set (int)  : template file choosen for lateral eye movement
%   eyemov_comp (int) : component for lateral eye movement in the template file
%
% Example function call:
% CORR_run(DATAPATH, studyname, blink_set, blink_comp, heart_set, heart_comp,eyemov_set, eyemov_comp)


% load STUDY
[STUDY ALLEEG] = pop_loadstudy('filename', [studyname, '.study'],'filepath', DATAPATH);
% setup parameters of the STUDY
CURRENTSTUDY = 1; 
EEG = ALLEEG; 
CURRENTSET = [1:length(EEG)];
[STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG);


% using CORRMAP algorithm to find correlation between ICA componenets and template 
% artefactual components in dataset
% 'th' : threshold for correlation
% 'ics' : independent component to work with
% 'pl' : plot the correlation map

% eye blinks 
[CORRMAP, STUDY, ALLEEG] = pop_corrmap(STUDY, ALLEEG, blink_set, blink_comp, 'th', '0.8', 'ics', 3, ...
    'pl', '2nd', 'title', 'plot', 'clname', '', 'badcomps', 'no', 'resetclusters', 'off');
% combining sets and ICs from second correlation map
eyeblink = [CORRMAP.output.sets{2} CORRMAP.output.ics{2}];
% variable to store the CORRMAP structure
Corr_eye = CORRMAP;

% heart beat 
[CORRMAP, STUDY, ALLEEG] = pop_corrmap(STUDY, ALLEEG, heart_set, heart_comp, 'th', '0.8', 'ics', 3, ...
    'pl', '2nd', 'title', 'plot', 'clname', '', 'badcomps', 'no', 'resetclusters', 'off');
% combining sets and ICs from second correlation map
heartbeat = [CORRMAP.output.sets{2} CORRMAP.output.ics{2}];
% variable to store the CORRMAP structure
Corr_heart = CORRMAP;

% lateral eye moveemnt 
[CORRMAP, STUDY, ALLEEG] = pop_corrmap(STUDY, ALLEEG, eyemov_set, eyemov_comp, 'th', '0.8', 'ics', 3, ...
    'pl', '2nd', 'title', 'plot', 'clname', '', 'badcomps', 'no', 'resetclusters', 'off');
% combining sets and ICs from second correlation map
eyemovement = [CORRMAP.output.sets{2} CORRMAP.output.ics{2}];
% variable to store the CORRMAP structure
Corr_lateye = CORRMAP;

% saving component info as MAT-file
save(fullfile(DATAPATH, 'components.mat'), 'eyeblink', 'heartbeat', 'eyemovement');
save(fullfile(DATAPATH, 'corrmap_info.mat'), 'Corr_eye', 'Corr_heart', 'Corr_lateye');