function [trialCount, corrEvent_left, corrEvent_right] = SATcorrect_trials(EEG)
% function [trialCount, corrEvent_left, corrEvent_right] = SATcorrect_trials(EEG)
%
% This function is used to retirieve the correct trial's event names from
% the EEG.events
% 
% Inputs:
%   EEG             (struct): EEGLab EEG data structure before any epoching or removing event names 
%
% Ouput:
%   trialCount      (int)   : total count of the trials
%   corrEvent_left  (cell)  : Event names of all the correct left trials
%   corrEvent_rigth (cell)  : Event names of all the correct right trials
%
% Example function call:
% [trialCount, corrEvent_left, corrEvent_right] = SATcorrect_trials(EEG)

% counter for trial
trialCount = 0;
% cell to store correct left events
corrEvent_left = {};
% cell to store correct right events
corrEvent_right = {};

% loop over events 
for idx = 1:numel(EEG.event)
    % check for response events 
    if strcmp(EEG.event(idx).type, 'response')
        % incrementing trial count 
        trialCount = trialCount + 1;
        
        % determining correct left trials
        % checking left trials
        if startsWith(EEG.event(idx-1).type, 'left')
            % checking for ascending trials
            if contains(EEG.event(idx-1).type, 'asc')
                % check if correct response
                if strcmp(EEG.event(idx+1).type, '1')
                    corrEvent_left = [corrEvent_left, EEG.event(idx-1).type];
                end 
            % checking for alternating trials
            elseif contains(EEG.event(idx-1).type, 'alt')
                % check if correct response
                if strcmp(EEG.event(idx+1).type, '2')
                    corrEvent_left = [corrEvent_left, EEG.event(idx-1).type];
                end 
            % checking for descending trials
            elseif contains(EEG.event(idx-1).type, 'dec')
                % check if correct response
                if strcmp(EEG.event(idx+1).type, '3')
                    corrEvent_left = [corrEvent_left, EEG.event(idx-1).type];
                end 
            end
            
        % determining correct right trials
        % checking right trials
        elseif startsWith(EEG.event(idx-1).type, 'right')
            % checking for ascending trials
            if contains(EEG.event(idx-1).type, 'asc')
                % check if correct response
                if strcmp(EEG.event(idx+1).type, '1')
                    corrEvent_right = [corrEvent_right, EEG.event(idx-1).type];
                end 
            % checking for alternating trials
            elseif contains(EEG.event(idx-1).type, 'alt')
                % check if correct response
                if strcmp(EEG.event(idx+1).type, '2')
                    corrEvent_right = [corrEvent_right, EEG.event(idx-1).type];
                end 
            % checking for descending trials
            elseif contains(EEG.event(idx-1).type, 'dec')
                % check if correct response
                if strcmp(EEG.event(idx+1).type, '3')
                    corrEvent_right = [corrEvent_right, EEG.event(idx-1).type];
                end 
            end  
        end 
    end 
end 





