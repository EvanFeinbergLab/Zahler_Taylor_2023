%% Load data and generate basic plots for SC-recipient PPRF opto-evoked gaze shifts 

BehaviorFilePath = '../Data/Figure 5 and supplements/SC-recipient PPRF stimulation';
Experiment_pprf = {
    
     %PPRF
    'DT689', [20211011]; %Example
    'DT724', [20211111]; 
    'DT725', [20211111];
    'DT728', [20211111];
    'SZ252', [20211116]; 
    'DT830', [20220321];
    'SZ350', [20220314];
    'SZ351', [20220314];
    'SZ352', [20220314];
    };

Variables = 'dlcPupilData,dlcLikelihoodData,dlcCamera2Data,dlcCamera2LikelihoodData,puffData,strain_gauge,config,experiment';
SessionNumber = [];

trial_param.saccade_threshold = 3;
trial_param.velocity_threshold = 100;
trial_param.frame_rate = 100;
trial_param.window_size = 101;
trial_param.align =51;
trial_param.response_window = trial_param.align-1:trial_param.align+10; 
trial_param.preAlignPeriod = 1:trial_param.align-2;
filter_param = struct();
filter_param.pupil_side = 2; % 0=left, 1=right, 2=mean
STRAIN_GAUGE_MEASUREMENT_POINT = 8;
head_mov_threshold = 0.2;

spont_head_mov_trial_count = [];
spont_total_trial_count = [];
for AnimalIdx = 1:9
    DateRange = [Experiment_pprf{AnimalIdx,2}];
    AnimalList = Experiment_pprf{AnimalIdx,1};
    BehaviorFiles=LoadBehavior('BehaviorFilePath', BehaviorFilePath, 'AnimalList', AnimalList, 'DateRange', DateRange,'Variables', Variables,'SessionNumber',SessionNumber);
    mtrx = combine_sessions(BehaviorFiles, trial_param, 1, 'FindSaccadesEye', 'mean');
    
    % Filter trials, plot opto-evoked traces, and compute saccade/head 
    % probabilities for summary statistics
    
    trial_param.starts = find(mtrx.optoTrial==1);
    [Filters,hankelMtrx] = Create_Filters(mtrx,trial_param,filter_param);
    FILTER = Filters.noNanFilter & ~Filters.saccadeBeforeAlignFilter & Filters.smoothStrainGaugeFilter & Filters.optoIntensityFilter==max(Filters.optoIntensityFilter); 
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    response_window_first_saccade = find_first_saccade_in_response_window(filterMtrx.saccade,trial_param.response_window);

    opto_saccade_amp = response_window_first_saccade;
    opto_init_pos = filterMtrx.initPupilPos;
    opto_end_pos = filterMtrx.finalPupilPos;
    opto_alignedPupil = filterMtrx.alignedPupil;
    opto_unalignedPupil = filterMtrx.unalignedPupil;
    opto_alignedStrainGauge = filterMtrx.alignedStrainGauge;
    opto_head_amp = opto_alignedStrainGauge(:,trial_param.align+STRAIN_GAUGE_MEASUREMENT_POINT); 

    prob_head_mov_pprf(AnimalIdx) = sum(abs(opto_head_amp)>head_mov_threshold)/length(opto_head_amp);
    prob_saccade_pprf(AnimalIdx) = sum(abs(opto_saccade_amp)>0)/length(opto_saccade_amp);
    prop_contra_saccade_pprf(AnimalIdx) = sum(opto_saccade_amp<0)/sum(abs(opto_saccade_amp)>0);
    
    if AnimalIdx == 1
        % Basic scatters and raw traces
        pupil_lim = [-25 25];
        head_lim = [-2 2];
        x1 = (-trial_param.align+1:trial_param.window_size-trial_param.align)/trial_param.frame_rate; % time series x-axis
        spont_color = [0.8 0.8 1];
        opto_color = [0 0 0];
        figure; set(gcf, 'Position',  [200, 100, 200, 200]);
        hold on; plot(x1,opto_unalignedPupil','color',opto_color);
        hold on; plot(x1,mean(opto_unalignedPupil),'r');
        title(AnimalList);
        xlabel('Time (s)'); ylabel('Pupil position (deg)');
        ylim(pupil_lim);
        figure; set(gcf, 'Position',  [200, 100, 200, 200]);
        hold on; plot(x1,opto_alignedStrainGauge','color',opto_color);
        hold on; plot(x1,mean(opto_alignedStrainGauge),'r');
        xlabel('Time (s)'); ylabel('Attempted head rot. (Z)');
        ylim(head_lim);
        figure; set(gcf, 'Position',  [200, 100, 200, 200]);
        hold on; scatter(opto_init_pos,opto_end_pos-opto_init_pos,10 ,opto_color,'filled');
        xlim(pupil_lim); ylim(pupil_lim);
        xlabel('Initial eye position (deg)'); ylabel('Saccade amplitude (deg)');
        figure; set(gcf, 'Position',  [200, 100, 200, 200]);
        hold on; scatter(opto_init_pos,opto_head_amp,10, opto_color,'filled');
        xlim(pupil_lim); ylim(head_lim);
        xlabel('Initial eye position (deg)'); ylabel('Head amplitude (Z)');
    end
    
    % Add additional filters to look at movements only. Use these data to
    % compute regression statistics
    
    %Saccades (PPRF only)
    FILTER = Filters.noNanFilter & Filters.saccadeAfterAlignFilter & ~Filters.saccadeBeforeAlignFilter & Filters.optoIntensityFilter==max(Filters.optoIntensityFilter); 
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    response_window_first_saccade = find_first_saccade_in_response_window(filterMtrx.saccade,trial_param.response_window);

    opto_saccade_amp = response_window_first_saccade;
    opto_init_pos = filterMtrx.initPupilPos;

    opto_saccade_amp_Zscore = (opto_saccade_amp-mean(opto_saccade_amp))/std(opto_saccade_amp);
    opto_init_pos_Zscore = (opto_init_pos-mean(opto_init_pos))/std(opto_init_pos);
    mdl = fitlm(opto_init_pos_Zscore,opto_saccade_amp_Zscore);
    slope_saccade_pprf(AnimalIdx) = table2array(mdl.Coefficients(2,1));
    R2_saccade_pprf(AnimalIdx) = mdl.Rsquared.Adjusted;

end



%% Load data and generate basic plots for SC-recipient Gi opto-evoked gaze shifts

BehaviorFilePath = '../Data/Figure 5 and supplements/SC-recipient Gi stimulation';
Experiment_gi = {
    
     %Gi
    'DT807', [20220224];
    'DT808', [20220224];    
    'DT809', [20220224];
    'DT811', [20220224]; %Example
    'DT819', [20220228];
    'DT820', [20220228];
    'DT821', [20220228];
    'DT825', [20220308];
    'DT826', [20220308];
    };

Variables = 'dlcPupilData,dlcLikelihoodData,dlcCamera2Data,dlcCamera2LikelihoodData,puffData,strain_gauge,config,experiment';
SessionNumber = [];

trial_param.saccade_threshold = 3;
trial_param.velocity_threshold = 100;
trial_param.frame_rate = 100;
trial_param.window_size = 101;
trial_param.align =51;
% If the saccade occurs within 1 frame of opto-onset, our saccade detection
% algorithm will mark onset as the frame before opto onset. Therefore, the
% response window must include the frame before opto onset
% (trial_param.align-1)
trial_param.response_window = trial_param.align-1:trial_param.align+10; 
trial_param.preAlignPeriod = 1:trial_param.align-2;
filter_param = struct();
filter_param.pupil_side = 2; % 0=left, 1=right, 2=mean
STRAIN_GAUGE_MEASUREMENT_POINT = 8;
head_mov_threshold = 0.2;

for AnimalIdx = 1:9
    DateRange = [Experiment_gi{AnimalIdx,2}];
    AnimalList = Experiment_gi{AnimalIdx,1};
    BehaviorFiles=LoadBehavior('BehaviorFilePath', BehaviorFilePath, 'AnimalList', AnimalList, 'DateRange', DateRange,'Variables', Variables,'SessionNumber',SessionNumber);
    mtrx = combine_sessions(BehaviorFiles, trial_param, 1, 'FindSaccadesEye', 'mean');
    
    % Filter trials, plot opto-evoked traces, and compute saccade/head 
    % probabilities for summary statistics

    trial_param.starts = find(mtrx.optoTrial==1);
    [Filters,hankelMtrx] = Create_Filters(mtrx,trial_param,filter_param);
    FILTER = Filters.noNanFilter & ~Filters.saccadeBeforeAlignFilter & Filters.smoothStrainGaugeFilter & Filters.optoIntensityFilter==max(Filters.optoIntensityFilter); 
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    response_window_first_saccade = find_first_saccade_in_response_window(filterMtrx.saccade,trial_param.response_window);

    opto_saccade_amp = response_window_first_saccade;
    opto_init_pos = filterMtrx.initPupilPos;
    opto_end_pos = filterMtrx.finalPupilPos;
    opto_alignedPupil = filterMtrx.alignedPupil;
    opto_unalignedPupil = filterMtrx.unalignedPupil;
    opto_alignedStrainGauge = filterMtrx.alignedStrainGauge;
    opto_head_amp = opto_alignedStrainGauge(:,trial_param.align+STRAIN_GAUGE_MEASUREMENT_POINT); 

    prob_head_mov_gi(AnimalIdx) = sum(abs(opto_head_amp)>head_mov_threshold)/length(opto_head_amp);
    prob_saccade_gi(AnimalIdx) = sum(abs(opto_saccade_amp)>0)/length(opto_saccade_amp);
    prop_contra_head_mov_gi(AnimalIdx) = sum(opto_head_amp<-head_mov_threshold)/sum(abs(opto_head_amp)>head_mov_threshold);

    if AnimalIdx == 4
        % Basic scatters and raw traces
        pupil_lim = [-25 25];
        head_lim = [-2 2];
        x1 = (-trial_param.align+1:trial_param.window_size-trial_param.align)/trial_param.frame_rate; % time series x-axis
        spont_color = [0.8 0.8 1];
        opto_color = [0 0 0];
        figure; set(gcf, 'Position',  [200, 100, 200, 200]);
        hold on; plot(x1,opto_unalignedPupil','color',opto_color);
        hold on; plot(x1,mean(opto_unalignedPupil),'r');
        title(AnimalList);
        xlabel('Time (s)'); ylabel('Pupil position (deg)');
        ylim(pupil_lim);
        figure; set(gcf, 'Position',  [200, 100, 200, 200]);
        hold on; plot(x1,opto_alignedStrainGauge','color',opto_color);
        hold on; plot(x1,mean(opto_alignedStrainGauge),'r');
        xlabel('Time (s)'); ylabel('Attempted head rot. (Z)');
        ylim(head_lim);
        figure; set(gcf, 'Position',  [200, 100, 200, 200]);
        hold on; scatter(opto_init_pos,opto_end_pos-opto_init_pos,10 ,opto_color,'filled');
        xlim(pupil_lim); ylim(pupil_lim);
        xlabel('Initial eye position (deg)'); ylabel('Saccade amplitude (deg)');
        figure; set(gcf, 'Position',  [200, 100, 200, 200]);
        hold on; scatter(opto_init_pos,opto_head_amp,10, opto_color,'filled');
        xlim(pupil_lim); ylim(head_lim);
        xlabel('Initial eye position (deg)'); ylabel('Head amplitude (Z)');
    end
    
    % Add additional filters to look at movements only. Use these data to
    % compute regression statistics
    
    %Head (Gi only)
    FILTER = Filters.noNanFilter & Filters.smoothStrainGaugeFilter & Filters.optoIntensityFilter==max(Filters.optoIntensityFilter); 
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);

    opto_alignedStrainGauge = filterMtrx.alignedStrainGauge;
    opto_head_amp = opto_alignedStrainGauge(:,trial_param.align+STRAIN_GAUGE_MEASUREMENT_POINT); 
    head_mov_FILTER = abs(opto_head_amp)>head_mov_threshold;
    opto_head_amp = opto_head_amp(head_mov_FILTER);
    opto_init_pos = opto_init_pos(head_mov_FILTER);

    opto_head_amp_Zscore = (opto_head_amp-mean(opto_head_amp))/std(opto_head_amp);
    opto_init_pos_Zscore = (opto_init_pos-mean(opto_init_pos))/std(opto_init_pos);
    mdl = fitlm(opto_init_pos_Zscore,opto_head_amp_Zscore);
    slope_head_gi(AnimalIdx) = table2array(mdl.Coefficients(2,1));
    R2_head_gi(AnimalIdx) = mdl.Rsquared.Adjusted;
    
end
%% Plot summary statistics

figure; set(gcf, 'Position',  [200, 100, 100, 200]);    
hold on; bar([mean(prob_head_mov_pprf),mean(prob_head_mov_gi)]);
hold on; scatter(ones(size(Experiment_pprf,1),1)+(rand(size(Experiment_pprf,1),1)*0.2-0.1), prob_head_mov_pprf,50,'k','filled')
hold on; scatter(ones(size(Experiment_gi,1),1).*2+(rand(size(Experiment_gi,1),1)*0.2-0.1), prob_head_mov_gi,50,'k','filled')
xlim([0 3]);
ylabel('P(head movement)');

figure; set(gcf, 'Position',  [200, 100, 100, 200]);
hold on; bar([mean(prob_saccade_pprf),mean(prob_saccade_gi)]);
hold on; scatter(ones(size(Experiment_pprf,1),1)+(rand(size(Experiment_pprf,1),1)*0.2-0.1), prob_saccade_pprf,50,'k','filled')
hold on; scatter(ones(size(Experiment_gi,1),1).*2+(rand(size(Experiment_gi,1),1)*0.2-0.1), prob_saccade_gi,50,'k','filled')
xlim([0 3]);
ylabel('P(saccade)');

figure; set(gcf, 'Position',  [200, 100, 100, 200]);
hold on; bar([mean(slope_saccade_pprf),mean(slope_head_gi)]);
hold on; scatter(ones(size(Experiment_pprf,1),1)+(rand(size(Experiment_pprf,1),1)*0.2-0.1), slope_saccade_pprf,50,'filled','k')
hold on; scatter(ones(size(Experiment_gi,1),1).*2+(rand(size(Experiment_gi,1),1)*0.2-0.1), slope_head_gi,50,'filled','k')
xlim([0 3]);
ylabel('slope');

figure; set(gcf, 'Position',  [200, 100, 100, 200]);
hold on; bar([mean(R2_saccade_pprf),mean(R2_head_gi)]);
hold on; scatter(ones(size(Experiment_pprf,1),1)+(rand(size(Experiment_pprf,1),1)*0.2-0.1), R2_saccade_pprf,50,'filled','k')
hold on; scatter(ones(size(Experiment_gi,1),1).*2+(rand(size(Experiment_gi,1),1)*0.2-0.1), R2_head_gi,50,'filled','k')
xlim([0 3]);
ylim([-0.1 0.6]);
ylabel('R2');


[~,p] = ttest2(prob_head_mov_gi,prob_head_mov_pprf,'Vartype', 'unequal')
[~,p] = ttest2(prob_saccade_gi,prob_saccade_pprf,'Vartype', 'unequal')
