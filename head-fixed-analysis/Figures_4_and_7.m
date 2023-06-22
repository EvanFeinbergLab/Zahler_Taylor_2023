%% Specify file paths and load data for ephys analyses


BehaviorFilePath = '../Data/Figure 4 and supplements/Tectoreticular ephys';
%BehaviorFilePath = '../Data/Figure 7 and supplements/PPRF ephys';
%BehaviorFilePath = '../Data/Figure 7 and supplements/Gi ephys';
%BehaviorFilePath = '../Data/Figure 7 and supplements/SC-recipient PPRF ephys';

Experiment = {
        % Tectoreticular animals/sessions
        'SZ204', [20211004]; % SC
        'SZ204', [20211005]; % SC
        'SZ205', [20211007]; % SC
        'SZ205', [20211008]; % SC
        'SZ205', [20211009]; % SC
        'SZ205', [20211011]; % SC
        'SZ215', [20211015]; % SC
        'SZ215', [20211017]; % SC
        'SZ253', [20211207]; % SC
%          
        % PPRF animals/sessions
%         'SZ261', [20220112]; % PPRF
%         'SZ302', [20220201]; % PPRF 
%         'SZ311', [20220201]; % PPRF 
%         'SZ308', [20220210]; % PPRF
%         'SZ362', [20220309]; % PPRF
%         'SZ361', [20220310]; % PPRF 
         
        % Gi animals/sessions
%         'SZ310', [20220209]; % Gi
%         'SZ338', [20220301]; % Gi
%         'SZ338', [20220302]; % Gi
%         'SZ345', [20220301]; % Gi
%         'SZ345', [20220302]; % Gi
%         'SZ362', [20220309]; % Gi

        % SC-recipient PPRF animals/session
%         'DTA90', [20230219]; 
%         'DTA92', [20230221]; 
%         'DTA92', [20230222]; 
%         'DTA94', [20230222]; 
%         'DTA95', [20230228]; 
%         'DTA95', [20230301]; 
%         'DTA95', [20230302]; 
%         'DTB21', [20230320]; 
%         'DTB26', [20230320]; 

    };

% variables to load
Variables = 'dlcPupilData,dlcLikelihoodData,dlcCamera2Data,dlcCamera2LikelihoodData,puffData,strain_gauge,config,experiment,sp,intan_times';


trial_param = struct();
trial_param.saccade_threshold = 3;
trial_param.velocity_threshold = 100;
trial_param.frame_rate = 100;
trial_param.window_size = 101;
trial_param.align = 51;
trial_param.response_window = trial_param.align:trial_param.align+10;
trial_param.preAlignPeriod = 1:trial_param.align-1;
trial_param.alignStrainGauge = trial_param.align;

filter_param = struct(); 
filter_param.smooth_strain_gauge_window = trial_param.align-10:trial_param.align;

STRAIN_GAUGE_MEASUREMENT_POINT = 8; % # of frames after airpuff to quantify head movement amplitude
head_mov_threshold = 0.2; % threshold for definining whether a head movement occured

mtrx_array = {};
sp_array = {};
intan_times_array = {};

for SessionIdx = 1:size(Experiment,1)
    DateRange = [Experiment{SessionIdx,2}];
    AnimalList = Experiment{SessionIdx,1};
    SessionNumber = []; %if empty, load all sessions on given date
    BehaviorFiles=LoadBehavior('BehaviorFilePath', BehaviorFilePath, 'AnimalList', AnimalList, 'DateRange', DateRange,'Variables', Variables,'SessionNumber',SessionNumber);
    mtrx_array{SessionIdx} = combine_sessions(BehaviorFiles, trial_param,1); % Convert behavior files into summary matrix
    sp_array{SessionIdx} = BehaviorFiles.sp; %ephys data
    intan_times_array{SessionIdx} = BehaviorFiles.intan_times; %digital inputs to intan board for synchronizing (airpuff onset, camera trigger onset)
    clear BehaviorFiles
end

%% Generate summary bar plots looking at sensory/motor tuning,: 

units = 2; % 0 = all units, 1 = untagged units, 2 = tagged units
period_of_interest = 1; % 0 = baseline, 1 = response, 2 = response-baseline
binSize = 0.001;
baseline_period = [-0.05 0];
response_period = [0 0.05];
baseline_length = abs(baseline_period(1)-baseline_period(2));
response_length = abs(response_period(1)-response_period(2));

p_saccade_amp = [];
p_head_amp = [];
p_end_pos = [];
p_puff = [];
slope_saccade_amp = [];
slope_head_amp = [];
slope_end_pos = [];
R2_saccade_amp = [];
R2_end_pos = [];
R2_head_amp = [];
puff_preference = [];
loc = [];

for SessionIdx = 1:size(Experiment,1)
    mtrx = mtrx_array{SessionIdx};
    sp = sp_array{SessionIdx};
    intan_times = intan_times_array{SessionIdx};
    
    tagged_cids = ID_tagged_Units(intan_times,sp,BehaviorFilePath,0);
    all_cids = sp.cids(sp.cgs==2);
    untagged_cids = setdiff(all_cids,tagged_cids);
    if units == 0 
        [spikeTimes, clu, cid] = filter_units(sp.st,sp.clu, all_cids);
    elseif units == 1
        [spikeTimes, clu, cid] = filter_units(sp.st,sp.clu, untagged_cids); 
    elseif units == 2
        [spikeTimes, clu, cid] = filter_units(sp.st,sp.clu, tagged_cids); 
    end
    trial_param.starts = find(mtrx.allPuff==1);
    [Filters,hankelMtrx] = Create_Filters(mtrx,trial_param,filter_param);
    
    % Look at all trials
    FILTER_puff = logical(Filters.includeAllFilter);  
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER_puff);
    puff_side = filterMtrx.rightPuff; %0 = left puff, 1 = right puff 
    
    % Look at saccade trials 
    FILTER_eye = logical(Filters.noNanFilter & Filters.saccadeAfterAlignFilter);        
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER_eye);
    response_window_first_saccade = find_first_saccade_in_response_window(filterMtrx.saccade,trial_param.response_window);
    saccade_amp = response_window_first_saccade;
    end_pos = filterMtrx.finalPupilPos;
        
    % Look at head movement trials (>0.2Z) not preceded by a head movement. 
    FILTER_head = logical(Filters.noNanFilter & Filters.smoothStrainGaugeFilter & abs(hankelMtrx.alignedStrainGauge(:,trial_param.align+STRAIN_GAUGE_MEASUREMENT_POINT))>head_mov_threshold);        
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER_head);
    alignedStrainGauge = filterMtrx.alignedStrainGauge;
    head_amp = alignedStrainGauge(:,trial_param.align+STRAIN_GAUGE_MEASUREMENT_POINT); 
    
    % Generate baseline and response psth. We will use this data to perform
    % statistical analyses for head and eye movements.    
    eventTimes = intan_times.allPuff(FILTER_puff);
    baseline_psth_puff = generate_psth(spikeTimes,clu,eventTimes,binSize,baseline_period);
    response_psth_puff = generate_psth(spikeTimes,clu,eventTimes,binSize,response_period);
    
    eventTimes = intan_times.allPuff(FILTER_eye);
    baseline_psth_eye = generate_psth(spikeTimes,clu,eventTimes,binSize,baseline_period);
    response_psth_eye = generate_psth(spikeTimes,clu,eventTimes,binSize,response_period);
    
    eventTimes = intan_times.allPuff(FILTER_head);
    baseline_psth_head = generate_psth(spikeTimes,clu,eventTimes,binSize,baseline_period);
    response_psth_head = generate_psth(spikeTimes,clu,eventTimes,binSize,response_period);
    
    % Perform statistics fo each neuron
    prev_length = length(p_saccade_amp);
    for neuron = 1:length(cid)
        response = sum(response_psth_puff(neuron,:,:),3)./response_length;
        baseline = sum(baseline_psth_puff(neuron,:,:),3)./baseline_length;
        if period_of_interest == 0
            data = baseline;
        elseif period_of_interest == 1
            data = response;
        elseif period_of_interest == 2
            data = response - baseline;
        end

        [~, p_puff(prev_length+neuron)] = ttest2(data(puff_side==0),data(puff_side==1));
        puff_preference(prev_length+neuron) = mean(data(puff_side==0))<mean(data(puff_side==1));
        
        response = sum(response_psth_eye(neuron,:,:),3)./response_length;
        baseline = sum(baseline_psth_eye(neuron,:,:),3)./baseline_length;
        if period_of_interest == 0
            data = baseline;
        elseif period_of_interest == 1
            data = response;
        elseif period_of_interest == 2
            data = response - baseline;
        end

        mdl_saccade_amp = fitlm(saccade_amp,data);  
        mdl_end_pos = fitlm(end_pos,data);  
        p_saccade_amp(prev_length+neuron) = table2array(mdl_saccade_amp.Coefficients(2,4));
        p_end_pos(prev_length+neuron) = table2array(mdl_end_pos.Coefficients(2,4));
        slope_saccade_amp(prev_length+neuron) = table2array(mdl_saccade_amp.Coefficients(2,1));
        slope_end_pos(prev_length+neuron) = table2array(mdl_end_pos.Coefficients(2,1));
        R2_saccade_amp(prev_length+neuron) = mdl_saccade_amp.Rsquared.Adjusted;
        R2_end_pos(prev_length+neuron) = mdl_end_pos.Rsquared.Adjusted;
        
        response = sum(response_psth_head(neuron,:,:),3)./response_length;
        baseline = sum(baseline_psth_head(neuron,:,:),3)./baseline_length;
        if period_of_interest == 0
            data = baseline;
        elseif period_of_interest == 1
            data = response;
        elseif period_of_interest == 2
            data = response - baseline;
        end

        mdl_head_amp = fitlm(head_amp,data);    
        p_head_amp(prev_length+neuron) = table2array(mdl_head_amp.Coefficients(2,4));
        slope_head_amp(prev_length+neuron) = table2array(mdl_head_amp.Coefficients(2,1));
        R2_head_amp(prev_length+neuron) = mdl_head_amp.Rsquared.Adjusted;

        if isfield(sp,'loc') 
            loc(prev_length+neuron,:) = sp.loc(logical(sp.cids==cid(neuron)),:);
        end
    end
end

p_val = 0.05; % signficance cut-off

total = length(p_saccade_amp);

contra_tuned_saccade_amp_bool = p_saccade_amp<p_val & slope_saccade_amp<0;
contra_tuned_saccade_amp = sum(contra_tuned_saccade_amp_bool);
ipsi_tuned_saccade_amp_bool = p_saccade_amp<p_val & slope_saccade_amp>0;
ipsi_tuned_saccade_amp = sum(ipsi_tuned_saccade_amp_bool);
tuned_saccade_amp_bool = p_saccade_amp<p_val;
not_tuned_saccade_amp = total-contra_tuned_saccade_amp-ipsi_tuned_saccade_amp;
binom_test_saccade_amp = myBinomTest(contra_tuned_saccade_amp,contra_tuned_saccade_amp+ipsi_tuned_saccade_amp,0.5,'two');

contra_tuned_head_amp_bool = p_head_amp<p_val & slope_head_amp<0;
contra_tuned_head_amp = sum(contra_tuned_head_amp_bool);
ipsi_tuned_head_amp_bool = p_head_amp<p_val & slope_head_amp>0;
ipsi_tuned_head_amp = sum(ipsi_tuned_head_amp_bool);
tuned_head_amp_bool = p_head_amp<p_val;
not_tuned_head_amp = total-contra_tuned_head_amp-ipsi_tuned_head_amp;
binom_test_head_amp = myBinomTest(contra_tuned_head_amp,contra_tuned_head_amp+ipsi_tuned_head_amp,0.5,'two');

contra_tuned_end_pos_bool = p_end_pos<p_val & slope_end_pos<0;
contra_tuned_end_pos = sum(contra_tuned_end_pos_bool);
ipsi_tuned_end_pos_bool = p_end_pos<p_val & slope_end_pos>0;
ipsi_tuned_end_pos = sum(ipsi_tuned_end_pos_bool);
tuned_end_pos_bool = p_end_pos<p_val;
not_tuned_end_pos = total-contra_tuned_end_pos-ipsi_tuned_end_pos;
binom_test_end_pos = myBinomTest(contra_tuned_end_pos,contra_tuned_end_pos+ipsi_tuned_end_pos,0.5,'two');

contra_tuned_puff_bool = puff_preference==0 & p_puff<p_val;
contra_tuned_puff = sum(contra_tuned_puff_bool);
ipsi_tuned_puff_bool = puff_preference==1 & p_puff<p_val;
ipsi_tuned_puff = sum(ipsi_tuned_puff_bool);
tuned_puff_bool = p_puff<p_val;
not_tuned_puff = total-contra_tuned_puff-ipsi_tuned_puff;
binom_test_puff = myBinomTest(contra_tuned_puff,contra_tuned_puff+ipsi_tuned_puff,0.5,'two');

% Generate bar graph
figure; set(gcf, 'Position',  [200, 100, 200, 200]);
bar([1;2;3;4], [contra_tuned_puff/total,ipsi_tuned_puff/total,not_tuned_puff/total;contra_tuned_saccade_amp/total,ipsi_tuned_saccade_amp/total,not_tuned_saccade_amp/total; contra_tuned_end_pos/total,ipsi_tuned_end_pos/total,not_tuned_end_pos/total; contra_tuned_head_amp/total,ipsi_tuned_head_amp/total,not_tuned_head_amp/total],'stacked'); 
xlim([0 5]);
xticklabels({});
fprintf('binomial tests\npuff side p = %d\nsaccade amp p = %d\nend pos p =%d\nhead amp p =%d\n',[binom_test_puff,binom_test_saccade_amp,binom_test_end_pos,binom_test_head_amp]);

fprintf('\nTuning counts\ncontra saccade:%d\nipsi saccade:%d\nnot tuned:%d\ncontra end pos:%d\nipsi end pos:%d\nnot tuned:%d\ncontra head:%d\nipsi head:%d\nnot tuned:%d\ncontra puff:%d\nipsi puff:%d\nnot tuned:%d\n',[contra_tuned_saccade_amp,ipsi_tuned_saccade_amp,not_tuned_saccade_amp,contra_tuned_end_pos,ipsi_tuned_end_pos,not_tuned_end_pos,contra_tuned_head_amp,ipsi_tuned_head_amp,not_tuned_head_amp,contra_tuned_puff,ipsi_tuned_puff,not_tuned_puff]);

% Generate matrix used for chi-square comparision between brain regions
chi_square_matrix_sc = [contra_tuned_saccade_amp,ipsi_tuned_saccade_amp,not_tuned_saccade_amp;contra_tuned_end_pos,ipsi_tuned_end_pos,not_tuned_end_pos;contra_tuned_head_amp,ipsi_tuned_head_amp,not_tuned_head_amp;contra_tuned_puff,ipsi_tuned_puff,not_tuned_puff];

%% Example unit PSTH and regression (aligned to puff onset)

period_of_interest = 1; % 0 = baseline, 1 = response, 2 = response-baseline
binSize = 0.001;
baseline_period = [-0.05 0];
response_period = [0 0.05];
baseline_length = abs(baseline_period(1)-baseline_period(2));
response_length = abs(response_period(1)-response_period(2));
psth_window = [-0.15 0.15];

% SC
% Example cell 1
example_session = 7;
example_neuron = 3018;

% Example cell 2
% example_session = 7;
% example_neuron = 2029;

% Example cell 3
% example_session = 8;
% example_neuron = 2039;


% PPRF
% Example cell 1
% example_session = 5;
% example_neuron = 3059;

% Example cell 2
% example_session = 5;
% example_neuron = 2011;

% Example cell 3
% example_session = 5;
% example_neuron = 3120;


% Gi
% Example cell 1
% example_session = 5;
% example_neuron = 1088;

% Example cell 2
% example_session = 5;
% example_neuron = 1112;

% Example cell 3
% example_session = 4;
% example_neuron = 4163;

% PPRF tagged
% Example cell 1
% example_session = 2;
% example_neuron = 2035;
% 
% Example cell 2
% example_session = 2;
% example_neuron = 2042;
% 
% Example cell 3
% example_session = 1;
% example_neuron = 1177;


mtrx = mtrx_array{example_session};
sp = sp_array{example_session};
intan_times = intan_times_array{example_session};
[spikeTimes, clu, cid] = filter_units(sp.st,sp.clu, example_neuron);

trial_param.starts = find(mtrx.allPuff==1);
[Filters,hankelMtrx] = Create_Filters(mtrx,trial_param,filter_param);

% Look at all trials
FILTER_puff = logical(Filters.includeAllFilter);  
filterMtrx = Create_filterMtrx(hankelMtrx,FILTER_puff);
puff_side = filterMtrx.rightPuff; %0 = left puff, 1 = right puff 

% Look at saccade trials not preceded by a saccade
FILTER_eye = logical(Filters.noNanFilter & Filters.saccadeAfterAlignFilter);        
filterMtrx = Create_filterMtrx(hankelMtrx,FILTER_eye);
response_window_first_saccade = find_first_saccade_in_response_window(filterMtrx.saccade,trial_param.response_window);
saccade_amp = response_window_first_saccade;
end_pos = filterMtrx.finalPupilPos;

% Look at head movement trials (>0.2Z) not preceded by a head movement. 
FILTER_head = logical(Filters.noNanFilter & Filters.smoothStrainGaugeFilter & abs(hankelMtrx.alignedStrainGauge(:,trial_param.align+STRAIN_GAUGE_MEASUREMENT_POINT))>head_mov_threshold);        
filterMtrx = Create_filterMtrx(hankelMtrx,FILTER_head);
alignedStrainGauge = filterMtrx.alignedStrainGauge;
head_amp = alignedStrainGauge(:,trial_param.align+STRAIN_GAUGE_MEASUREMENT_POINT); 

%Generate baseline and response psth. We will use this data to perform
%regression analyses for head and eye movements.
eventTimes = intan_times.allPuff(FILTER_eye);
baseline_psth_eye = generate_psth(spikeTimes,clu,eventTimes,binSize,baseline_period);
response_psth_eye = generate_psth(spikeTimes,clu,eventTimes,binSize,response_period);

eventTimes = intan_times.allPuff(FILTER_head);
baseline_psth_head = generate_psth(spikeTimes,clu,eventTimes,binSize,baseline_period);
response_psth_head = generate_psth(spikeTimes,clu,eventTimes,binSize,response_period);

%Regression for eye movements
response = sum(response_psth_eye(1,:,:),3)./response_length;
baseline = sum(baseline_psth_eye(1,:,:),3)./baseline_length;
if period_of_interest == 0
    data = baseline;
elseif period_of_interest == 1
    data = response;
elseif period_of_interest == 2
    data = response - baseline;
end

mdl_saccade_amp = fitlm(saccade_amp,data);  
mdl_end_pos = fitlm(end_pos,data);  
p_saccade_amp = table2array(mdl_saccade_amp.Coefficients(2,4));
p_end_pos = table2array(mdl_end_pos.Coefficients(2,4));
slope_saccade_amp = table2array(mdl_saccade_amp.Coefficients(2,1));
slope_end_pos = table2array(mdl_end_pos.Coefficients(2,1));
R2_saccade_amp = mdl_saccade_amp.Rsquared.Adjusted;
R2_end_pos = mdl_end_pos.Rsquared.Adjusted;

figure; set(gcf, 'Position',  [200, 100, 200, 200]); 
scatter(saccade_amp,data,'filled');
coords = predict(mdl_saccade_amp,[-12; 12]);
hold on; plot([-12; 12], coords);
ylim([0, 200]);
xlim([-12, 12])
fprintf('\nSaccade amp:\nSlope = %d\np = %d\nR2 = %d\n',[slope_saccade_amp,p_saccade_amp,R2_saccade_amp]);

figure; set(gcf, 'Position',  [200, 100, 200, 200]); 
scatter(end_pos,data,'filled');
coords = predict(mdl_end_pos,[-12; 12]);
hold on; plot([-12; 12], coords);
ylim([0, 200]);
xlim([-12, 12])
fprintf('\nSaccade end pos:\nSlope = %d\np = %d\nR2 = %d\n',[slope_end_pos,p_end_pos,R2_end_pos]);

%Regression for head movements
response = sum(response_psth_head(1,:,:),3)./response_length;
baseline = sum(baseline_psth_head(1,:,:),3)./baseline_length;
if period_of_interest == 0
    data = baseline;
elseif period_of_interest == 1
    data = response;
elseif period_of_interest == 2
    data = response - baseline;
end

mdl_head_amp = fitlm(head_amp,data);    
p_head_amp = table2array(mdl_head_amp.Coefficients(2,4));
slope_head_amp = table2array(mdl_head_amp.Coefficients(2,1));
R2_head_amp = mdl_head_amp.Rsquared.Adjusted;

figure; set(gcf, 'Position',  [200, 100, 200, 200]);
scatter(head_amp,data,'filled');
coords = predict(mdl_head_amp,[-1; 1]);
hold on; plot([-1; 1], coords);
ylim([0, 200]);
xlim([-1, 1])
fprintf('\nHead amp:\nSlope = %d\np = %d\nR2 = %d\n',[slope_head_amp,p_head_amp,R2_head_amp]);

% Bin data so that we can visualize PSTHs for puffs and movements based on 
% left vs right. 
puff_side_bins = puff_side+1;
num_bins=2;
saccade_amp_bins = createBins(saccade_amp,[-30, 30],num_bins);
num_bins=2;
saccade_end_bins = createBins(end_pos,[-30 30],num_bins);
num_bins=2;
head_amp_bins = createBins(head_amp,[-30, 30],num_bins);

% Set limit for PSTHs. Need to adjust this depending on the unit.
ylimits = [0 250];

eventTimes = intan_times.allPuff(FILTER_puff);
trial_idx = [];
groups = [];
for i = 1:2
    groups = [groups; ones(sum(puff_side_bins==i),1) * i];
    trial_idx = [trial_idx; find(puff_side_bins==i)];
end
groups(groups==0) = [];
psthViewer(spikeTimes, clu, eventTimes(trial_idx), psth_window,groups);
subplot(3,1,1); ylim(ylimits)

eventTimes = intan_times.allPuff(FILTER_eye);
trial_idx = [];
groups = [];
for i = 1:2
    groups = [groups; ones(sum(saccade_amp_bins==i),1) * i];
    trial_idx = [trial_idx; find(saccade_amp_bins==i)];
end
groups(groups==0) = [];
psthViewer(spikeTimes, clu, eventTimes(trial_idx), psth_window,groups);
subplot(3,1,1); ylim(ylimits)

eventTimes = intan_times.allPuff(FILTER_eye);
trial_idx = [];
groups = [];
for i = 1:2
    groups = [groups; ones(sum(saccade_end_bins==i),1) * i];
    trial_idx = [trial_idx; find(saccade_end_bins==i)];
end
groups(groups==0) = [];
psthViewer(spikeTimes, clu, eventTimes(trial_idx), psth_window,groups);
subplot(3,1,1); ylim(ylimits)

eventTimes = intan_times.allPuff(FILTER_head);
trial_idx = [];
groups = [];
for i = 1:2
    groups = [groups; ones(sum(head_amp_bins==i),1) * i];
    trial_idx = [trial_idx; find(head_amp_bins==i)];
end
groups(groups==0) = [];
psthViewer(spikeTimes, clu, eventTimes(trial_idx), psth_window,groups);
subplot(3,1,1); ylim(ylimits)
