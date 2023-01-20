counter=0;
muscleChannelsAll=[];

for iSN=1:length(Subject);
    load([DUKEDIR '/' Subject(iSN).Name '/muscleChannels.mat'])
    muscleChannels=muscleChannels+counter;
    muscleChannelsAll=cat(2,muscleChannelsAll,muscleChannels);
    counter=counter+length(Subject(iSN).goodChannels);
end