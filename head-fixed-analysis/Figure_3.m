%% Load data and generate basic plots for Tectoreticular opto-evoked gaze shifts 

BehaviorFilePath = '../Data/Figure 3 and supplements/Tectoreticular terminal stimulation';
Experiment_pmrf = {
    'SZ050', [20210509]; %Example
    'SZ052', [20210509]; 
    'SZ388', [20220420]; 
    'SZ390', [20220420]; 
    'SZ394', [20220420]; 
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

for AnimalIdx = 1:5
    DateRange = [Experiment_pmrf{AnimalIdx,2}];
    AnimalList = Experiment_pmrf{AnimalIdx,1};
    BehaviorFiles=LoadBehavior('BehaviorFilePath', BehaviorFilePath, 'AnimalList', AnimalList, 'DateRange', DateRange,'Variables', Variables,'SessionNumber',SessionNumber);
    mtrx = combine_sessions(BehaviorFiles, trial_param, 1, 'FindSaccadesEye', 'mean');
        
    trial_param.starts = find(mtrx.optoTrial==1);
    [Filters,hankelMtrx] = Create_Filters(mtrx,trial_param,filter_param);

    %Saccades
    FILTER = Filters.noNanFilter & Filters.saccadeAfterAlignFilter & ~Filters.saccadeBeforeAlignFilter & Filters.optoIntensityFilter==max(Filters.optoIntensityFilter); 
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    response_window_first_saccade = find_first_saccade_in_response_window(filterMtrx.saccade,trial_param.response_window);

    opto_saccade_amp = response_window_first_saccade;
    opto_init_pos_saccade_trials = filterMtrx.initPupilPos;
    opto_end_pos = filterMtrx.finalPupilPos;
    opto_unalignedPupil = filterMtrx.unalignedPupil;

    opto_saccade_amp_Zscore = (opto_saccade_amp-mean(opto_saccade_amp))/std(opto_saccade_amp);
    opto_init_pos_Zscore = (opto_init_pos_saccade_trials-mean(opto_init_pos_saccade_trials))/std(opto_init_pos_saccade_trials);   
    mdl = fitlm(opto_init_pos_Zscore,opto_saccade_amp_Zscore);
    slope_saccade_pmrf(AnimalIdx) = table2array(mdl.Coefficients(2,1));
    R2_saccade_pmrf(AnimalIdx) = mdl.Rsquared.Adjusted;

    %Head
    FILTER = Filters.noNanFilter & Filters.smoothStrainGaugeFilter & Filters.optoIntensityFilter==max(Filters.optoIntensityFilter) & abs(hankelMtrx.alignedStrainGauge(:,trial_param.align+STRAIN_GAUGE_MEASUREMENT_POINT))>head_mov_threshold; 
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    
    opto_init_pos_head_trials = filterMtrx.initPupilPos;
    opto_alignedStrainGauge = filterMtrx.alignedStrainGauge;
    opto_head_amp = opto_alignedStrainGauge(:,trial_param.align+STRAIN_GAUGE_MEASUREMENT_POINT); 
    head_mov_FILTER = abs(opto_head_amp)>head_mov_threshold;
    opto_head_amp = opto_head_amp(head_mov_FILTER);
    opto_init_pos_head_trials = opto_init_pos_head_trials(head_mov_FILTER);

    opto_head_amp_Zscore = (opto_head_amp-mean(opto_head_amp))/std(opto_head_amp);
    opto_init_pos_Zscore = (opto_init_pos_head_trials-mean(opto_init_pos_head_trials))/std(opto_init_pos_head_trials);
    mdl = fitlm(opto_init_pos_Zscore,opto_head_amp_Zscore);
    slope_head_mov_pmrf(AnimalIdx) = table2array(mdl.Coefficients(2,1));
    R2_head_mov_pmrf(AnimalIdx) = mdl.Rsquared.Adjusted;

    
    if AnimalIdx == 1
        % Basic scatters and raw traces
        pupil_lim = [-25 25];
        head_lim = [-4 4];
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
        hold on; scatter(opto_init_pos_saccade_trials,opto_saccade_amp,10 ,opto_color,'filled');
        xlim(pupil_lim); ylim(pupil_lim);
        xlabel('Initial eye position (deg)'); ylabel('Saccade amplitude (deg)');
        figure; set(gcf, 'Position',  [200, 100, 200, 200]);
        hold on; scatter(opto_init_pos_head_trials,opto_head_amp,10, opto_color,'filled');
        xlim(pupil_lim); ylim(head_lim);
        xlabel('Initial eye position (deg)'); ylabel('Head amplitude (Z)');
   end

end

%Plot summary data
figure; set(gcf, 'Position',  [200, 100, 100, 200]);
hold on; plot([slope_saccade_pmrf;slope_head_mov_pmrf],'k');
hold on; scatter(ones(size(Experiment_pmrf,1),1)+(rand(size(Experiment_pmrf,1),1)*0.2-0.1),slope_saccade_pmrf ,50,'filled','k')
hold on; scatter(ones(size(Experiment_pmrf,1),1).*2+(rand(size(Experiment_pmrf,1),1)*0.2-0.1),slope_head_mov_pmrf  ,50,'filled','k')
xlim([0 3]);
ylim([-1 0.4]);
ylabel('slope');
xticklabels({'','Eyes','Head',''});

figure; set(gcf, 'Position',  [200, 100, 100, 200]);
hold on; plot([R2_saccade_pmrf;R2_head_mov_pmrf],'k');
hold on; scatter(ones(size(Experiment_pmrf,1),1)+(rand(size(Experiment_pmrf,1),1)*0.2-0.1),R2_saccade_pmrf ,50,'filled','k')
hold on; scatter(ones(size(Experiment_pmrf,1),1).*2+(rand(size(Experiment_pmrf,1),1)*0.2-0.1),R2_head_mov_pmrf ,50,'filled','k')
xlim([0 3]);
ylim([-0.1 1]);
ylabel('R2');
xticklabels({'','Eyes','Head',''});

fprintf('\nSaccade statistics:\nMean slope: %d\nStd dev slope: %d\nMean R2: %d\nStd dev R2: %d',[mean(slope_saccade_pmrf),std(slope_saccade_pmrf),mean(R2_saccade_pmrf),std(R2_saccade_pmrf)]);
fprintf('\nHead statistics:\nMean slope: %d\nStd dev slope: %d\nMean R2: %d\nStd dev R2: %d\nMean probability: %d\nStd probability: %d\n',[mean(slope_head_mov_pmrf),std(slope_head_mov_pmrf),mean(R2_head_mov_pmrf),std(R2_head_mov_pmrf)]);
[~, p]= ttest(slope_saccade_pmrf,slope_head_mov_pmrf)
[~, p]= ttest(R2_saccade_pmrf,R2_head_mov_pmrf)

%% Load data and generate basic plots for tectoreticular Halo inhibition

BehaviorFilePath = '../Data/Figure 3 and supplements/Tectoreticular inhibition (hsv-cre+halo)/Halo';

Experiment = {  
        'SZ186', [20210830:20210924];
        'SZ187', [20210830:20210924];
        'SZ229', [20211026:20211129];  
        'SZ230', [20211026:20211129]; 
        'SZ231', [20211026:20211129]; 
        'SZ271', [20211213:20211221];
        'SZ273', [20211213:20211221];
        'SZ275', [20211213:20211221];
        'SZ276', [20211213:20211221];
        'SZ277', [20211213:20211221];
        'SZ278', [20211213:20211221];
        'SZ283', [20211213:20211221];
    };

Variables = 'dlcPupilData,dlcLikelihoodData,dlcCamera2Data,dlcCamera2LikelihoodData,puffData,strain_gauge,config,experiment';
SessionNumber = [];

trial_param.saccade_threshold = 3;
trial_param.velocity_threshold = 100;
trial_param.frame_rate = 100;
trial_param.window_size = 101;
trial_param.align =51;
trial_param.response_window = trial_param.align:trial_param.align+10; 
trial_param.preAlignPeriod = 1:trial_param.align-1;
filter_param = struct();
filter_param.pupil_side = 2; % 0=left, 1=right, 2=mean
STRAIN_GAUGE_MEASUREMENT_POINT = 8;
head_mov_threshold = 0.2;

for AnimalIdx = 1:12
    DateRange = [Experiment{AnimalIdx,2}];
    AnimalList = Experiment{AnimalIdx,1};
    BehaviorFiles=LoadBehavior('BehaviorFilePath', BehaviorFilePath, 'AnimalList', AnimalList, 'DateRange', DateRange,'Variables', Variables,'SessionNumber',SessionNumber);
    mtrx = combine_sessions(BehaviorFiles, trial_param, 1, 'FindSaccadesEye', 'mean');
        
    trial_param.starts = find(mtrx.allPuff==1);
    [Filters,hankelMtrx] = Create_Filters(mtrx,trial_param,filter_param);
    
    %SACCADES
    
    % Left control trials 
    FILTER = Filters.noNanFilter & Filters.leftPuffFilter & ~Filters.optoTrialFilter & Filters.saccadeAfterAlignFilter; 
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    response_window_first_saccade = find_first_saccade_in_response_window(filterMtrx.saccade,trial_param.response_window);
    left_control_saccade_amp = response_window_first_saccade;
    left_control_end_pos = filterMtrx.finalPupilPos;
    % Left opto trials
    FILTER = Filters.noNanFilter & Filters.leftPuffFilter & Filters.optoTrialFilter & Filters.saccadeAfterAlignFilter; 
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    response_window_first_saccade = find_first_saccade_in_response_window(filterMtrx.saccade,trial_param.response_window);
    left_opto_saccade_amp = response_window_first_saccade;
    left_opto_end_pos = filterMtrx.finalPupilPos;
    % Right control trials 
    FILTER = Filters.noNanFilter & Filters.rightPuffFilter & ~Filters.optoTrialFilter & Filters.saccadeAfterAlignFilter; 
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    response_window_first_saccade = find_first_saccade_in_response_window(filterMtrx.saccade,trial_param.response_window);
    right_control_saccade_amp = response_window_first_saccade;
    right_control_end_pos = filterMtrx.finalPupilPos;
    % Right opto trials
    FILTER = Filters.noNanFilter & Filters.rightPuffFilter & Filters.optoTrialFilter & Filters.saccadeAfterAlignFilter; 
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    response_window_first_saccade = find_first_saccade_in_response_window(filterMtrx.saccade,trial_param.response_window);
    right_opto_saccade_amp = response_window_first_saccade;
    right_opto_end_pos = filterMtrx.finalPupilPos;   
    
    mean_left_control_end_pos(AnimalIdx) = mean(left_control_end_pos);
    mean_left_opto_end_pos(AnimalIdx) = mean(left_opto_end_pos);
    mean_left_control_saccade_amp(AnimalIdx) = mean(left_control_saccade_amp);
    mean_left_opto_saccade_amp(AnimalIdx) = mean(left_opto_saccade_amp);
    mean_right_control_end_pos(AnimalIdx) = mean(right_control_end_pos);
    mean_right_opto_end_pos(AnimalIdx) = mean(right_opto_end_pos);
    mean_right_control_saccade_amp(AnimalIdx) = mean(right_control_saccade_amp);
    mean_right_opto_saccade_amp(AnimalIdx) = mean(right_opto_saccade_amp);
    
    %HEAD MOVEMENTS
    
    % Left control trials 
    FILTER = Filters.noNanFilter & Filters.leftPuffFilter & ~Filters.optoTrialFilter; 
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    tmp_head_amp = filterMtrx.alignedStrainGauge(:,trial_param.align+STRAIN_GAUGE_MEASUREMENT_POINT); 
    left_control_head_amp = tmp_head_amp(abs(tmp_head_amp)>head_mov_threshold);
    % Left opto trials 
    FILTER = Filters.noNanFilter & Filters.leftPuffFilter & Filters.optoTrialFilter; 
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    tmp_head_amp = filterMtrx.alignedStrainGauge(:,trial_param.align+STRAIN_GAUGE_MEASUREMENT_POINT); 
    left_opto_head_amp = tmp_head_amp(abs(tmp_head_amp)>head_mov_threshold);
    % Right control trials 
    FILTER = Filters.noNanFilter & Filters.rightPuffFilter & ~Filters.optoTrialFilter; 
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    tmp_head_amp = filterMtrx.alignedStrainGauge(:,trial_param.align+STRAIN_GAUGE_MEASUREMENT_POINT); 
    right_control_head_amp = tmp_head_amp(abs(tmp_head_amp)>head_mov_threshold);
    % Right opto trials 
    FILTER = Filters.noNanFilter & Filters.rightPuffFilter & Filters.optoTrialFilter; 
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    tmp_head_amp = filterMtrx.alignedStrainGauge(:,trial_param.align+STRAIN_GAUGE_MEASUREMENT_POINT); 
    right_opto_head_amp = tmp_head_amp(abs(tmp_head_amp)>head_mov_threshold);
    
    mean_left_control_head_amp(AnimalIdx) = mean(left_control_head_amp);
    mean_left_opto_head_amp(AnimalIdx) = mean(left_opto_head_amp);
    mean_right_control_head_amp(AnimalIdx) = mean(right_control_head_amp);
    mean_right_opto_head_amp(AnimalIdx) = mean(right_opto_head_amp);  
end

figure; set(gcf, 'Position',  [200, 100, 100, 200]);
hold on; bar(mean(mean_left_control_end_pos-mean_left_opto_end_pos));
hold on; scatter(ones(size(Experiment,1),1)+(rand(size(Experiment,1),1)*0.2-0.1),mean_left_control_end_pos-mean_left_opto_end_pos,50,'filled','k')
xlim([0 3]);
ylim([-3 3]);

figure; set(gcf, 'Position',  [200, 100, 100, 200]);
hold on; bar(mean(mean_left_control_head_amp-mean_left_opto_head_amp));
hold on; scatter(ones(size(Experiment,1),1)+(rand(size(Experiment,1),1)*0.2-0.1),mean_left_control_head_amp-mean_left_opto_head_amp,50,'filled','k')
xlim([0 3]);
ylim([-0.5 0.5])

fprintf('\n left whisker airpuff endpoint (CONTROL) = %g +- %g std \n', [mean(mean_left_control_end_pos), std(mean_left_control_end_pos)])
fprintf('\n left whisker airpuff endpoint (OPTO) = %g +- %g std \n', [mean(mean_left_opto_end_pos), std(mean_left_opto_end_pos)])
[~, p] = ttest(mean_left_control_end_pos,mean_left_opto_end_pos);
fprintf('\n left whisker airpuff endpoint CONTROL vs. OPTO p = %g \n', [p]);

fprintf('\n left whisker airpuff head amp (CONTROL) = %g +- %g std \n', [mean(mean_left_control_head_amp), std(mean_left_control_head_amp)])
fprintf('\n left whisker airpuff head amp (OPTO) = %g +- %g std \n', [mean(mean_left_opto_head_amp), std(mean_left_opto_head_amp)])
[~, p] = ttest(mean_left_control_head_amp,mean_left_opto_head_amp);
fprintf('\n left whisker airpuff head amp CONTROL vs. OPTO p = %g \n', [p]);
