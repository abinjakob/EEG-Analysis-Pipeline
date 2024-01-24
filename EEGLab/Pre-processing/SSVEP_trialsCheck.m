L20count = 0;
L15count = 0;
R20count = 0;
R15count = 0;

for c = 1:size(EEG.event,2)
    if strcmp(EEG.event(c).type, events{1})
        L20count = L20count + 1;
    elseif strcmp(EEG.event(c).type, events{2})
        L15count = L15count + 1;
    elseif strcmp(EEG.event(c).type, events{3})
        R20count = R20count + 1;
    elseif strcmp(EEG.event(c).type, events{4})
        R15count = R15count + 1;
    end
end 

display(['Left 20 = ' num2str(L20count)]);
display(['Left 15 = ' num2str(L15count)]);
display(['Right 20 = ' num2str(R20count)]);
display(['Right 15 = ' num2str(R15count)]);
display(['Total Trials = ' num2str(L20count + L15count + R20count + R15count)]);