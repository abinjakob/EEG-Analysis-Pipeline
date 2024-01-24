event_sel = 1;
stims = {'Left 20Hz', 'Left 15Hz', 'Right 20Hz', 'Right 15Hz'};
% extracting the event data 
EEG_temp = pop_selectevent(EEG, 'type', events{event_sel},'renametype', events{event_sel}, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
% parameters for pwelch
window_length = EEG_temp.pnts;
overlap = window_length / 2;

% calculating power 
% loop over channels 
for iChan = 1:size(EEG_temp.data,1)
    % loop over trials 
    for iTrial = 1:size(EEG_temp.data,3)
        % computing psd usign pwelch
        [pxx, f] = pwelch(EEG_temp.data(iChan,:,iTrial), hamming(window_length), overlap, 2^nextpow2(window_length*4), EEG_temp.srate);
        % computing psd usign pwelch
        pxx_all(:,iChan,iTrial) = pxx;
    end
end 

% plotting single trial power
chan = [5 39 49 50];
% figure(3), clf
% % loop over trials
% h = plot(f, squeeze(pxx_all(:,chan,:)));
% hold on
% plot(f, mean(pxx_all(:,chan,:),3),'k', 'linew', 1.2)
% set(h, 'Color', [200,200,200]/255)
% set(gca, 'xlim', [5 50], 'ylim', [0 10]);
% xlabel('frequency [Hz]')
% ylabel('power [a.u.]')
% title(['PSD of channel ' EEG.chanlocs(chan).labels ' in Attend ' stims{event_sel} ' condition'])

% plotting 
figure(4), clf
% imagesc(f, [], squeeze(pxx_all(:,chan,:))')
imagesc(f, [], squeeze(mean(pxx_all(:,chan,:),2))')
set(gca, 'xlim', [5 45], 'clim', [0 20])
xlabel('frequency [Hz]')
ylabel('trials')
% title(['PSD of channel ' EEG.chanlocs(chan).labels ' in Attend ' stims{event_sel} ' condition'])
title(['PSD of Occipital Channels in Attend ' stims{event_sel} ' condition'])
colorbar