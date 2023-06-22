function [psth] = generate_psth(spikeTimes,clu,eventTimes,binSize,window)
%generate_psth generates a psth aligned to the specified event
%   spikeTimes      timestamps of all spikes
%   clu             unit ID for each spike
%   eventTimes      timestamps of events you want your psth aligned to
%   binSize         size of each bin in ms (usually 1 ms)
%   window          period surrounding each event you'd like extracted

clusterIDs = sort(unique(clu));
binBorders = window(1):binSize:window(2);
numBins = length(binBorders)-1;
psth = zeros(length(clusterIDs),length(eventTimes),numBins);

for clusterIndex = 1:length(clusterIDs)
    st = spikeTimes(clu==clusterIDs(clusterIndex));

    for r = 1:length(eventTimes)
        [n, ~] = histdiff(st, eventTimes(r), binBorders);
        psth(clusterIndex,r,:) = n;
    end
end

end
