function [spikeTimes_filtered,clu_filtered,clusterID_filtered] = filter_units(spikeTimes,clu,filter_units)
%Inputs
    %spikeTimes: Spike timstamps
    %clu: spike cluster id
    %filter_idx: cluster id of units you want to keep
    
%Outputs
    %spikeTime_filtered: spiketimes of units you want to keep
    %clu_filtered: corresponding cluster ids of units you want to keep
    spikeTimes_filtered = [];
    clu_filtered = [];
    for i = 1:length(filter_units)
        spikeTimes_filtered = [spikeTimes_filtered; spikeTimes(clu == filter_units(i))];
        clu_filtered = [clu_filtered; clu(clu == filter_units(i))];
    end
    [spikeTimes_filtered, sort_idx] = sort(spikeTimes_filtered);
    clu_filtered = clu_filtered(sort_idx);
    clusterID_filtered = sort(unique(clu_filtered));
    
end

