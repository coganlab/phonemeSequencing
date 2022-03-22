function [ieegPowerSeries,ieegFilt,ieegPower] = extractPower(ieegSplit,fs,Frange,time,tw)
    for i = 1 : size(ieegSplit,1)
        ieegFiltChannel = eegfilt(squeeze(ieegSplit(i,:,:)),fs,Frange(1),Frange(2),0,200);
        ieegFilt(i,:,:) = ieegFiltChannel;
        for tr = 1:size(ieegSplit,2)
            %ieegPower(i,tr) = mean(log10(ieegFiltChannel(tr,time>=tw(1)&time<=tw(2)).^2));
            ieegPower(i,tr) = mean((ieegFiltChannel(tr,time>=tw(1)&time<=tw(2))));
            ieegPowerSeries(i,tr,:) = (abs(hilbert(ieegFiltChannel(tr,:))));
        end
    end    
end