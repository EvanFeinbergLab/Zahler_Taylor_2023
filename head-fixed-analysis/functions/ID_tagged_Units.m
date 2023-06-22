function [ag_cid] = ID_tagged_Units(intan_times,sp,BehaviorFilePath,plot_me)
%ID_tagged_units Identify units that are responsive to opto stim. 
%   intan_times: contains timestamps for LED onset
%   sp: contains spike times for each unit
%   BehaviorFilePath: root of ephys folder. ex: *_shank1_ephys
%   plot_me: if 1, plot psth and waveforms for each optotagged unit
    try
        %My stimulation protocol involves 20 pulses at session start and 20
        %pulses at session end
        eventTimes = intan_times.LED_on([1:20, end-19:end]); 

        %Filter "good" units identified during spike sorting
        [g_spikeTimes, g_clu, g_cid] = filter_units(sp.st,sp.clu,sp.cids(sp.cgs==2));

        binSize = 0.001;
        window = [-0.01 0.01];
        psth = generate_psth(g_spikeTimes,g_clu,eventTimes,binSize,window); 

        p_increase = zeros(size(psth,1),1);
        response_rate = zeros(size(psth,1),1);
        jitter = zeros(size(psth,1),1);
        for neuron = 1:size(psth,1)
            wn = window * 1000;   % convert window size to ms
            dt = binSize * 1000;   % convert time resolution of bin raster to ms

            baseline = sum(psth(neuron,:,1:10),3); % sum baseline spikes for each trial
            response = sum(psth(neuron,:,11:20),3); % sum response period spikes for each trial
            [~, p_increase(neuron)] = ttest(baseline,response,'tail','left'); %ttest looking for increase compared to baseline

            %Find spike latency for each trial
            first_spike = [];
            for trial = 1:size(psth,2)
                trial_spikes = find(psth(neuron,trial,11:15)>0);
                if ~isempty(trial_spikes)
                    first_spike(trial) = trial_spikes(1);
                else
                    first_spike(trial) = 0; % no response
                end
            end
            
            response_rate(neuron) = 1-sum(first_spike==0)/size(psth,2); % subtract fraction of no-response trials from 1 to get response rate
            jitter(neuron) = std(first_spike(first_spike~=0)); 
        end

        % metrics for determining opto-responsiveness. 
        antidromic = p_increase<0.001 & response_rate>0.8 & jitter<2;
        [ag_spikeTimes, ag_clu, ag_cid] = filter_units(g_spikeTimes,g_clu,g_cid(antidromic));
        
       

        if plot_me == 1
            %View antidromically activated neurons
            
            % PSTH of opto-evoked spikes for tagged neurons
            window = [-0.02 0.02];
            binSize = 0.001;
            psth = generate_psth(ag_spikeTimes,ag_clu,eventTimes,binSize,window); 
            x_idx = linspace(window(1),window(2),(window(2)-window(1))/binSize);
            shank_folder_list = floor(ag_cid/1000);
            unit_number_list = rem(ag_cid,1000);
            
            for neuron = 1:length(ag_cid)
                figure; set(gcf, 'Position',  [200, 100, 100, 200]); 
                title(sprintf('shank %d unit %d',shank_folder_list(neuron), unit_number_list(neuron)))
                for trial = 1 : size(psth,2)
                    x1 = x_idx((psth(neuron,trial,:)~=0));
                    y1 = ones(length(find(psth(neuron,trial,:)~=0)),1).*trial; 
                    hold on; scatter(x1,y1,5,'filled','k')
                end
                xlim([x_idx(1),x_idx(end)]);
                ylim([1 size(psth,2)]);
            end
            
            try
                %Compute correlation between mean waveform and opto-evoked waveform
                for neuron = 1:length(ag_cid)
                    shank_folder = sprintf([BehaviorFilePath,'/',BehaviorFilePath(end-16:end),'_shank%d_ephys'], shank_folder_list(neuron));
                    unit_number = double(unit_number_list(neuron));
                    sp_shank = loadKSdir(shank_folder);
                    gwfparams. dataDir = shank_folder;
                    gwfparams.fileName = 'temp_wh.dat';
                    gwfparams.dataType = 'int16';
                    gwfparams.nCh = 32;
                    gwfparams.wfWin = [-20 20];
                    gwfparams.nWf = 200;
                    gwfparams.spikeTimes = round(sp_shank.st(sp_shank.clu==unit_number).*30000);
                    gwfparams.spikeClusters = unit_number.*ones(size(find(sp_shank.clu==unit_number)));

                    wf = getWaveForms(gwfparams);

                    [~, max_wf_idx] = max(range(wf.waveFormsMean(1,:,:),3));

                    if plot_me == 1
                        figure; set(gcf, 'Position',  [200, 100, 300, 1000]);
                        for i = 1:32
                            subtightplot(16,2,i)
                            control_wf = squeeze(wf.waveForms(1,:,i,:))';
                            hold on; plot(control_wf,'k');
                        end
                    end
                    control_wf_max = squeeze(wf.waveFormsMean(1,max_wf_idx(1),:))';

                    spike_idx = round(sp_shank.st(sp_shank.clu==unit_number).*30000);
                    opto_LED_times = intan_times.LED_on([1:20, 21:40])*30000;
                    opto_spikes = [];
                    for i = 1:length(opto_LED_times)
                        tmp = spike_idx(spike_idx>opto_LED_times(i) & spike_idx<=opto_LED_times(i)+150);
                        opto_spikes = [opto_spikes; tmp];
                    end

                    gwfparams.spikeTimes = opto_spikes;
                    gwfparams.nWf = length(opto_spikes);
                    gwfparams.spikeClusters = unit_number.*ones(size(opto_spikes));
                    wf = getWaveForms(gwfparams);

                    if plot_me == 1
                        for i = 1:32
                            subtightplot(16,2,i)
                            opto_wf = squeeze(wf.waveFormsMean(:,i,:))';
                            hold on; plot(opto_wf,'color','r');
                            ylim([-3000 2000])
                            xticklabels([]);
                            yticklabels([]);
                        end
                        subtightplot(16,2,1)
                        title(sprintf('shank %d unit %d',shank_folder_list(neuron), unit_number_list(neuron)))
                    end
                    opto_wf_max = squeeze(wf.waveFormsMean(1,max_wf_idx(1),:))';

                    mdl = fitlm(control_wf_max(15:30),opto_wf_max(15:30));
                    fprintf('shank %d unit %d: corr. coef. %d',shank_folder_list(neuron), unit_number_list(neuron),mdl.Rsquared.Adjusted);  
                end
            catch
            end
            
        end
        
        
        
    catch
        fprintf('No antidromic units');
    end
  
end

