%% Load data

BehaviorFilePath = '../Data/Figure 2 and supplements/SC Electrical Stim';
Experiment = {
    'DT731', [20211108:20211108], [01];
    'DT751', [20211124:20211124], [01];
    'DT752', [20211124:20211124], [01];
    'DT759', [20211216:20211216], [01];
    'DT780', [20220121:20220121], [01];
    'DT798', [20220128:20220128], [01]; % example animal
    'DT799', [20220129:20220129], [01];
    'DT803', [20220215:20220215], [01];
    'DT822', [20220222:20220222], [03]; 
    }; %#ok<*NBRAK2> 

% Set trial parameters for event detection (these are used by combine_sessions_freely_moving)
trial_param.saccade_threshold = 3; % deg
trial_param.velocity_threshold = 100; % deg/sec
trial_param.frame_rate = 100; % Hz
mtrx_array = {};

Variables = 'dlcPupilData,dlcLikelihoodData,dlcCamera2Data,dlcCamera2LikelihoodData,puffData,strain_gauge,config,experiment';  % variables to load

% Load data and create finalized time series
for AnimalIdx = 1:size(Experiment,1)
    DateRange = [Experiment{AnimalIdx,2}];
    AnimalList = Experiment{AnimalIdx,1};
    SessionNumber = Experiment{AnimalIdx,3};
    BehaviorFiles=LoadBehavior('BehaviorFilePath', BehaviorFilePath, 'AnimalList', AnimalList, 'DateRange', DateRange,'Variables', Variables,'SessionNumber',SessionNumber);
    mtrx_array{AnimalIdx} = combine_sessions(BehaviorFiles, trial_param, 1, 'FindSaccadesEye', 'mean'); %#ok<*SAGROW> 
end

% Define colors
c_pupil = [127, 63, 152]/255;
c_head = [247, 148, 29]/255;
c_pupil_light = (1 - c_pupil)*0.5 + c_pupil;
c_head_light = (1 - c_head)*0.5 + c_head;

%% CALCULATE TRIAL VARIABLES, PLOT EXAMPLE ANIMALS

example_animal_idx = 6;
% example_animal_idx = 1:9;

trial_param.window_size = 101;
trial_param.align =51;
trial_param.response_window = trial_param.align-1:trial_param.align+10;
trial_param.preAlignPeriod = 1:trial_param.align-2;
trial_param.alignStrainGauge = trial_param.align;
trial_param.window_xaxis = ((1:trial_param.window_size) - trial_param.align)*1/trial_param.frame_rate;
filter_param = struct();
filter_param.pupil_side = 2; % 0=left, 1=right, 2=mean
STRAIN_GAUGE_MEASUREMENT_POINT = 8;

endpoint_mean = [];
saccade_amp_mean = [];
endpoint_std = [];
saccade_amp_std = [];
saccade_frac_left = [];

for AnimalIdx = 1:size(Experiment,1)
    mtrx = mtrx_array{AnimalIdx};
    
    shank_filters = {'Filters.leftPuff2Filter', 'Filters.leftPuffFilter', 'Filters.rightPuffFilter', 'Filters.rightPuff2Filter'}; % leftPuff2=most anterior site, rightPuff2=most posterior site
    num_shanks = numel(shank_filters);

    for SHANK = 1:num_shanks
        trial_param.starts = find(mtrx.allPuff==1);
        [Filters,hankelMtrx] = Create_Filters(mtrx,trial_param,filter_param);
        FILTER = Filters.noNanFilter & Filters.smoothStrainGaugeFilter & eval(shank_filters{SHANK}) & Filters.saccadeAfterAlignFilter;
        filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
        shankFilterMtrx = filterMtrx;

        shank(SHANK).init_pos = shankFilterMtrx.unalignedPupil(:, trial_param.align);
        shank(SHANK).saccade_amp = find_first_saccade_in_response_window(shankFilterMtrx.saccade, trial_param.response_window);
        shank(SHANK).end_pos = shank(SHANK).init_pos + shank(SHANK).saccade_amp;
        shank(SHANK).unalignedPupil = shankFilterMtrx.unalignedPupil;
        shank(SHANK).alignedStrainGauge = shankFilterMtrx.alignedStrainGauge;
        shank(SHANK).head_amp = shankFilterMtrx.alignedStrainGauge(:, trial_param.align + STRAIN_GAUGE_MEASUREMENT_POINT);
        
        % Means and variance
        endpoint_mean(AnimalIdx,SHANK) = mean(shank(SHANK).end_pos);
        endpoint_std(AnimalIdx,SHANK) = std(shank(SHANK).end_pos);
        saccade_amp_mean(AnimalIdx,SHANK) = mean(shank(SHANK).saccade_amp);
        saccade_amp_std(AnimalIdx,SHANK) = std(shank(SHANK).saccade_amp);
        head_amp_mean(AnimalIdx,SHANK) = mean(shank(SHANK).head_amp);
        head_amp_std(AnimalIdx,SHANK) = std(shank(SHANK).head_amp);
        
        % movement direction proportions
        saccade_frac_left(AnimalIdx,SHANK) = sum(shank(SHANK).saccade_amp<0)/numel(shank(SHANK).saccade_amp);
        head_frac_left(AnimalIdx,SHANK) = sum(shank(SHANK).head_amp<0)/numel(shank(SHANK).head_amp);

        % Regressions
        % standardize coefficients with 'normalize'
        mdl = fitlm( normalize(shank(SHANK).init_pos), normalize(shank(SHANK).head_amp) );
        pupil_v_head_r2(AnimalIdx,SHANK) = mdl.Rsquared.Ordinary;
        pupil_v_head_coef(AnimalIdx,SHANK) = mdl.Coefficients{2,1};

        mdl = fitlm( normalize(shank(SHANK).init_pos), normalize(shank(SHANK).saccade_amp) );
        pupil_v_pupil_r2(AnimalIdx,SHANK) = mdl.Rsquared.Ordinary;
        pupil_v_pupil_coef(AnimalIdx,SHANK) = mdl.Coefficients{2,1};

        mdl = fitlm( normalize(shank(SHANK).init_pos), normalize(shank(SHANK).end_pos) );
        pupil_v_pupil_endpoint_r2(AnimalIdx,SHANK) = mdl.Rsquared.Ordinary;
        pupil_v_pupil_endpoint_coef(AnimalIdx,SHANK) = mdl.Coefficients{2,1};


        % SACCADE PROBABILITIES
        % get counts for all trials and bin centers
        num_bins = 4;
        pupil_bin_range = [-10, 10];
        pupil_bin_edges = linspace(pupil_bin_range(1), pupil_bin_range(2), num_bins+1);

        FILTER = Filters.noNanFilter & Filters.smoothStrainGaugeFilter & eval(shank_filters{SHANK});
        filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
        [N_pupil_all, pupil_bin_edges] = histcounts(filterMtrx.unalignedPupil(:,trial_param.align), pupil_bin_edges);
        pupil_bin_width = diff(pupil_bin_edges);
        pupil_bin_width = pupil_bin_width(1);
        pupil_bin_centers = pupil_bin_edges(1:end-1)+pupil_bin_width/2;
        shank(SHANK).probabilities.saccade_bin_centers(AnimalIdx, :) = pupil_bin_centers;
        
        % all saccades
        FILTER = Filters.noNanFilter & Filters.saccadeAfterAlignFilter & Filters.smoothStrainGaugeFilter & eval(shank_filters{SHANK});
        filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
        N_saccade = histcounts(filterMtrx.unalignedPupil(:,trial_param.align), pupil_bin_edges); % after filtering, see how many trials fall into initial eye position bins
        shank(SHANK).probabilities.saccade(AnimalIdx, :) = N_saccade./N_pupil_all;
    
        % left saccades
        FILTER = Filters.noNanFilter & Filters.leftSaccadeAfterAlignFilter & Filters.smoothStrainGaugeFilter & eval(shank_filters{SHANK});
        filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
        N_saccade = histcounts(filterMtrx.unalignedPupil(:,trial_param.align), pupil_bin_edges); % after filtering, see how many trials fall into initial eye position bins
        shank(SHANK).probabilities.saccade_left(AnimalIdx, :) = N_saccade./N_pupil_all;
    
        % right saccades
        FILTER = Filters.noNanFilter & Filters.rightSaccadeAfterAlignFilter & Filters.smoothStrainGaugeFilter & eval(shank_filters{SHANK});
        filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
        N_saccade = histcounts(filterMtrx.unalignedPupil(:,trial_param.align), pupil_bin_edges); % after filtering, see how many trials fall into initial eye position bins
        shank(SHANK).probabilities.saccade_right(AnimalIdx, :) = N_saccade./N_pupil_all;


        % HEAD MOVEMENT PROBABILITIES
        % get counts for all trials and bin centers
        num_bins = 4;
        pupil_bin_range = [-10, 10];
        pupil_bin_edges = linspace(pupil_bin_range(1), pupil_bin_range(2), num_bins+1);

        FILTER = Filters.noNanFilter & Filters.smoothStrainGaugeFilter & eval(shank_filters{SHANK});
        filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
        [N_all, pupil_bin_edges] = histcounts(filterMtrx.unalignedPupil(:,trial_param.align), pupil_bin_edges);
        pupil_bin_width = diff(pupil_bin_edges);
        pupil_bin_width = pupil_bin_width(1);
        pupil_bin_centers = pupil_bin_edges(1:end-1)+pupil_bin_width/2;
        shank(SHANK).probabilities.head_mov_bin_centers(AnimalIdx, :) = pupil_bin_centers;
        
        % all head movements
        FILTER = Filters.noNanFilter & abs(hankelMtrx.alignedStrainGauge(:,trial_param.align + STRAIN_GAUGE_MEASUREMENT_POINT)) > 0.2 & Filters.smoothStrainGaugeFilter & eval(shank_filters{SHANK});
        filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
        N_saccade = histcounts(filterMtrx.unalignedPupil(:,trial_param.align), pupil_bin_edges); % after filtering, see how many trials fall into initial eye position bins
        shank(SHANK).probabilities.head_mov(AnimalIdx, :) = N_saccade./N_all;
    
        % left head movements
        FILTER = Filters.noNanFilter & hankelMtrx.alignedStrainGauge(:,trial_param.align + STRAIN_GAUGE_MEASUREMENT_POINT) < -0.2 & Filters.smoothStrainGaugeFilter & eval(shank_filters{SHANK});
        filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
        N_saccade = histcounts(filterMtrx.unalignedPupil(:,trial_param.align), pupil_bin_edges); % after filtering, see how many trials fall into initial eye position bins
        shank(SHANK).probabilities.head_mov_left(AnimalIdx, :) = N_saccade./N_all;
    
        % right head movements
        FILTER = Filters.noNanFilter & hankelMtrx.alignedStrainGauge(:,trial_param.align + STRAIN_GAUGE_MEASUREMENT_POINT) > 0.2 & Filters.smoothStrainGaugeFilter & eval(shank_filters{SHANK});
        filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
        N_saccade = histcounts(filterMtrx.unalignedPupil(:,trial_param.align), pupil_bin_edges); % after filtering, see how many trials fall into initial eye position bins
        shank(SHANK).probabilities.head_mov_right(AnimalIdx, :) = N_saccade./N_all;

    end
    
    % Attempted head movement traces
    if any(AnimalIdx == example_animal_idx)
        headlimits = [-6 6];
        figure; set(gcf, 'Position',  [200, 300, 800, 200]); hold on
        for SHANK = 1:num_shanks
            subplot(1,4,SHANK);
            plot(trial_param.window_xaxis, shank(SHANK).alignedStrainGauge, '-', 'Color', c_head);
            ylim(headlimits);
            ylabel('Head force. (Z)')
            axis square
        end
    end

    % Pupil traces
    if any(AnimalIdx == example_animal_idx)
        pupillimits = [-30 30];
        figure; set(gcf, 'Position',  [200, 300, 800, 200]); hold on
        for SHANK = 1:num_shanks
            subplot(1,4,SHANK);
            plot(trial_param.window_xaxis, shank(SHANK).unalignedPupil, '-', 'Color', c_pupil);
            ylim(pupillimits);
            ylabel('Pupil position (°)')
            axis square
        end
    end

    % Initial pupil angle vs. head displacement
    if any(AnimalIdx == example_animal_idx)
        figure; set(gcf, 'Position',  [200, 300, 800, 200]); hold on
        for SHANK = 1:num_shanks
            subplot(1,4,SHANK); hold on
            lim_size = 30;
            scatter(shank(SHANK).init_pos, shank(SHANK).head_amp, 20, c_head)
            xlim([-lim_size, lim_size])
            ylim([-4, 4])
            vline(0, 'k--')
            hline(0, 'k--')
            text_spacing = 4/4;
            mdl = fitlm(shank(SHANK).init_pos, shank(SHANK).head_amp);
            plot([-lim_size, lim_size], predict(mdl,[-4; 4]), 'k-')
            text(-lim_size+text_spacing,6-text_spacing*1, sprintf('R2 = %0.2f', mdl.Rsquared.Ordinary))
            text(-lim_size+text_spacing,6-text_spacing*2, sprintf('Slope = %0.2f', mdl.Coefficients{2,1}))
            xlabel('Initial pupil angle (°)')
            ylabel('Head disp. (°)')
            axis square
        end
    end

    % Initial pupil angle vs. saccade displacement
    if any(AnimalIdx == example_animal_idx)
        figure; set(gcf, 'Position',  [200, 300, 800, 200]); hold on
        for SHANK = 1:num_shanks
            subplot(1,4,SHANK); hold on
            lim_size = 30;
            scatter(shank(SHANK).init_pos, shank(SHANK).saccade_amp, 20, c_pupil)
            xlim([-lim_size, lim_size])
            ylim([-lim_size, lim_size])
            vline(0, 'k--')
            hline(0, 'k--')
            text_spacing = lim_size/6;
            mdl = fitlm(shank(SHANK).init_pos, shank(SHANK).saccade_amp);
            plot([-lim_size, lim_size], predict(mdl,[-lim_size; lim_size]), 'k-')
            text(-lim_size+text_spacing,lim_size-text_spacing*1, sprintf('R2 = %0.2f', mdl.Rsquared.Ordinary))
            text(-lim_size+text_spacing,lim_size-text_spacing*2, sprintf('Slope = %0.2f', mdl.Coefficients{2,1}))
            xlabel('Initial pupil angle (°)')
            ylabel('Saccade disp. (°)')
            axis square
        end
    end


end

%% PLOT PUPIL MEAN ENDPOINTS AND AMPLITUDES

example_animal_idx = 6;

% Plot mean pupil endpoints by shank
figure; set(gcf, 'Position',  [200, 300, 200, 200]); hold on
plot(1:4, endpoint_mean', '-', 'Color', c_pupil_light);  % plot all traces
plot(1:4, mean(endpoint_mean, 1)', '-', 'LineWidth', 3, 'Color', c_pupil) % plot the average trace
ylim([-15 15])
xlim([0.5 4.5])
xticks(1:4)
xlabel('Shank')
ylabel('Saccade endpoint (°)')

% Plot mean pupil displacements by shank
figure; set(gcf, 'Position',  [200, 300, 200, 200]); hold on
plot(1:4, saccade_amp_mean', '-', 'Color', c_pupil_light);  % plot all traces
plot(1:4, mean(saccade_amp_mean, 1)', '-', 'LineWidth', 3, 'Color', c_pupil) % plot the average trace
ylim([-15 15])
xlim([0.5 4.5])
xticks(1:4)
xlabel('Shank')
ylabel('Saccade displacement (°)')

% Plot mean head displacement by shank
figure; set(gcf, 'Position',  [200, 300, 200, 200]); hold on
plot(1:4, head_amp_mean', '-', 'Color', c_head_light); % plot all traces
plot(1:4, mean(head_amp_mean, 1)', '-', 'LineWidth', 3, 'Color', c_head) % plot the average trace
ylim([-4, 4])
xlim([0.5 4.5])
xticks(1:4)
xlabel('Shank')
ylabel('Attempted head disp (Z)')


% Calculate statistics
tmp_mean = mean(endpoint_mean, 1);
tmp_std = std(endpoint_mean, 1);
fprintf('\nMean saccade endpoint:\n')
fprintf('Shank 1: %f +- %f  (mean +- SD)\n', tmp_mean(1), tmp_std(1))
fprintf('Shank 2: %f +- %f  (mean +- SD)\n', tmp_mean(2), tmp_std(2))
fprintf('Shank 3: %f +- %f  (mean +- SD)\n', tmp_mean(3), tmp_std(3))
fprintf('Shank 4: %f +- %f  (mean +- SD)\n', tmp_mean(4), tmp_std(4))

tmp_mean = mean(saccade_amp_mean, 1);
tmp_std = std(saccade_amp_mean, 1);
fprintf('\nMean saccade displacement:\n')
fprintf('Shank 1: %f +- %f (mean +- SD)\n', tmp_mean(1), tmp_std(1))
fprintf('Shank 2: %f +- %f (mean +- SD)\n', tmp_mean(2), tmp_std(2))
fprintf('Shank 3: %f +- %f (mean +- SD)\n', tmp_mean(3), tmp_std(3))
fprintf('Shank 4: %f +- %f (mean +- SD)\n', tmp_mean(4), tmp_std(4))
fprintf('Overall: %f +- %f (mean +- SD)\n', mean(saccade_amp_mean(:)), std(saccade_amp_mean(:)))

tmp_mean = mean(head_amp_mean, 1);
tmp_std = std(head_amp_mean, 1);
fprintf('\nMean head displacement:\n')
fprintf('Shank 1: %f +- %f (mean +- SD)\n', tmp_mean(1), tmp_std(1))
fprintf('Shank 2: %f +- %f (mean +- SD)\n', tmp_mean(2), tmp_std(2))
fprintf('Shank 3: %f +- %f (mean +- SD)\n', tmp_mean(3), tmp_std(3))
fprintf('Shank 4: %f +- %f (mean +- SD)\n', tmp_mean(4), tmp_std(4))
fprintf('Overall: %f +- %f (mean +- SD)\n', mean(head_amp_mean(:)), std(head_amp_mean(:)))


tmp_x = repmat(1:4, size(head_amp_mean,1), 1);
tmp_x = tmp_x(:);
tmp_y = head_amp_mean(:);
mdl = fitlm(tmp_x, tmp_y);
fprintf('\nRegression of shank against mean head displacement:\n')
fprintf('R2: %f\n', mdl.Rsquared.Ordinary)
fprintf('p: %f\n', mdl.Coefficients{"x1", "pValue"})

tmp_x = repmat(1:4, size(saccade_amp_mean,1), 1);
tmp_x = tmp_x(:);
tmp_y = saccade_amp_mean(:);
mdl = fitlm(tmp_x, tmp_y);
fprintf('\nRegression of shank against mean saccade displacement:\n')
fprintf('R2: %f\n', mdl.Rsquared.Ordinary)
fprintf('p: %f\n', mdl.Coefficients{"x1", "pValue"})

tmp_x = repmat(1:4, size(endpoint_mean,1), 1);
tmp_x = tmp_x(:);
tmp_y = endpoint_mean(:);
mdl = fitlm(tmp_x, tmp_y);
fprintf('\nRegression of shank against mean saccade endpoint:\n')
fprintf('R2: %f\n', mdl.Rsquared.Ordinary)
fprintf('p: %f\n', mdl.Coefficients{"x1", "pValue"})


%% PLOT STANDARD DEVIATIONS OF PUPIL ENDPOINTS AND DISPLACEMENTS

% Plot overall saccade displacement/endpoint variability
figure; set(gcf, 'Position',  [200, 300, 200, 200]); hold on
bar(1, mean(mean(saccade_amp_std, 2)), 'FaceColor', c_pupil_light)
bar(2, mean(mean(endpoint_std, 2)), 'FaceColor', c_pupil_light)
plot(1:2, [mean(saccade_amp_std, 2), mean(endpoint_std, 2)], 'o-', 'Color', [0 0 0], 'MarkerSize',4);
xticks([1, 2])
xticklabels({'Disp.', 'Endpoint'})
xtickangle(20)
ylabel('Variability (SD)')

% STATISTICS
fprintf('\nSaccade endpoint variability (SD) (all shanks): %f +- %f (mean +- SD)\n', mean(mean(endpoint_std, 2)), std(mean(endpoint_std, 2)))
fprintf('Saccade amplitude variability (SD) (all shanks): %f +- %f (mean +- SD)\n', mean(mean(saccade_amp_std, 2)), std(mean(saccade_amp_std, 2)))
[~, p] = ttest2(mean(endpoint_std, 2), mean(saccade_amp_std, 2));
fprintf('t-test: Saccade amplitude variability vs saccade endpoint variability p: %f\n', p)


%% PLOT REGRESSION STATS FOR HEAD/SACCADE AMP AGAINST INITIAL PUPIL POSITION

% Plot R2 and slope for regression of head displacement against initial
% pupil position
figure; set(gcf, 'Position',  [200, 300, 200, 200]); hold on
bar(1:4, mean(pupil_v_head_r2, 1), 'FaceColor', c_head_light)
plot(1:4, pupil_v_head_r2, 'ko-', 'MarkerSize',4);
ylim([0 1])
ylabel('R2')
xticks(1:4)
xlabel('Shank')

figure; set(gcf, 'Position',  [200, 300, 200, 200]); hold on
bar(1:4, mean(pupil_v_head_coef, 1), 'FaceColor', c_head_light)
plot(1:4, pupil_v_head_coef, 'ko-', 'MarkerSize',4);
ylim([-1.1 0.2])
ylabel('Slope')
xticks(1:4)
xlabel('Shank')


% Plot R2 and slope for regression of saccade displacement against initial pupil position
figure; set(gcf, 'Position',  [200, 300, 200, 200]); hold on
bar(1:4, mean(pupil_v_pupil_r2, 1), 'FaceColor', c_pupil_light)
plot(1:4, pupil_v_pupil_r2, 'ko-', 'MarkerSize',4);
ylim([0 1])
ylabel('R2')
xticks(1:4)
xlabel('Shank')

figure; set(gcf, 'Position',  [200, 300, 200, 200]); hold on
bar(1:4, mean(pupil_v_pupil_coef, 1), 'FaceColor', c_pupil_light)
plot(1:4, pupil_v_pupil_coef, 'ko-', 'MarkerSize',4);
ylim([-1.1 0.2])
ylabel('Slope')
xticks(1:4)
xlabel('Shank')


% Plot R2 values together (with bounded line)
figure; set(gcf, 'Position',  [200, 300, 200, 200]); hold on
boundedline(1:4, mean(pupil_v_pupil_r2, 1), std(pupil_v_pupil_r2, 1),'cmap', c_pupil);
boundedline(1:4, mean(pupil_v_head_r2, 1), std(pupil_v_head_r2, 1),'cmap', c_head);
ylim([0 1])
ylabel('Variance explained by initial pupil pos (R2)')
xticks(1:4)
xlabel('Shank')

% Plot slope values together (with bounded line)
figure; set(gcf, 'Position',  [200, 300, 200, 200]); hold on
boundedline(1:4, mean(pupil_v_pupil_coef, 1), std(pupil_v_pupil_coef, 1),'cmap', c_pupil);
boundedline(1:4, mean(pupil_v_head_coef, 1), std(pupil_v_head_coef, 1),'cmap', c_head);
% ylim([0 1])
ylabel('Slope')
xticks(1:4)
xlabel('Shank')


% STATISTICS

tmp_mean = mean(pupil_v_head_r2, 1);
tmp_std = std(pupil_v_head_r2, 1);
fprintf('\nInitial eye vs head amp (R2):\n')
fprintf('Shank 1: %f +- %f (mean +- SD)\n', tmp_mean(1), tmp_std(1))
fprintf('Shank 2: %f +- %f (mean +- SD)\n', tmp_mean(2), tmp_std(2))
fprintf('Shank 3: %f +- %f (mean +- SD)\n', tmp_mean(3), tmp_std(3))
fprintf('Shank 4: %f +- %f (mean +- SD)\n', tmp_mean(4), tmp_std(4))
fprintf('Overall: %f +- %f (mean +- SD)\n', mean(pupil_v_head_r2(:)), std(pupil_v_head_r2(:)))

tmp_mean = mean(pupil_v_pupil_r2, 1);
tmp_std = std(pupil_v_pupil_r2, 1);
fprintf('\nInitial eye vs saccade amp (R2):\n')
fprintf('Shank 1: %f +- %f (mean +- SD)\n', tmp_mean(1), tmp_std(1))
fprintf('Shank 2: %f +- %f (mean +- SD)\n', tmp_mean(2), tmp_std(2))
fprintf('Shank 3: %f +- %f (mean +- SD)\n', tmp_mean(3), tmp_std(3))
fprintf('Shank 4: %f +- %f (mean +- SD)\n', tmp_mean(4), tmp_std(4))
fprintf('Overall: %f +- %f (mean +- SD)\n', mean(pupil_v_pupil_r2(:)), std(pupil_v_pupil_r2(:)))


%% PLOT FRACTION OF CONTRA SACCADES AND HEAD MOVEMENTS

% SACCADES
figure; set(gcf, 'Position',  [200, 300, 200, 200]); hold on
bar(1:4, mean(saccade_frac_left, 1), 'FaceColor', c_pupil_light)
plot(1:4, saccade_frac_left, 'ko-', 'MarkerSize', 4);
ylim([0 1])
ylabel('Frac. contra saccades')
xticks(1:4)
xlabel('Shank')

% ATTEMPTED HEAD
figure; set(gcf, 'Position',  [200, 300, 200, 200]); hold on
bar(1:4, mean(head_frac_left, 1), 'FaceColor', c_head_light)
plot(1:4, head_frac_left, 'ko-', 'MarkerSize', 4);
ylim([0 1])
ylabel('Frac. contra head mov.')
xticks(1:4)
xlabel('Shank')


%STATISTICS
tmp_mean = mean(saccade_frac_left, 1);
tmp_std = std(saccade_frac_left, 1);
fprintf('\nFraction leftward (contraversive) pupil movements:\n')
fprintf('Shank 1: %f +- %f (mean +- SD)\n', tmp_mean(1), tmp_std(1))
fprintf('Shank 2: %f +- %f (mean +- SD)\n', tmp_mean(2), tmp_std(2))
fprintf('Shank 3: %f +- %f (mean +- SD)\n', tmp_mean(3), tmp_std(3))
fprintf('Shank 4: %f +- %f (mean +- SD)\n', tmp_mean(4), tmp_std(4))
fprintf('Overall: %f +- %f (mean +- SD)\n', mean(saccade_frac_left(:)), std(saccade_frac_left(:)))

tmp_mean = mean(head_frac_left, 1);
tmp_std = std(head_frac_left, 1);
fprintf('\nFraction leftward (contraversive) head movements:\n')
fprintf('Shank 1: %f +- %f (mean +- SD)\n', tmp_mean(1), tmp_std(1))
fprintf('Shank 2: %f +- %f (mean +- SD)\n', tmp_mean(2), tmp_std(2))
fprintf('Shank 3: %f +- %f (mean +- SD)\n', tmp_mean(3), tmp_std(3))
fprintf('Shank 4: %f +- %f (mean +- SD)\n', tmp_mean(4), tmp_std(4))
fprintf('Overall: %f +- %f (mean +- SD)\n', mean(head_frac_left(:)), std(head_frac_left(:)))


%% Predict shank position with multinomial logistic regression

animals = 1:8; % omitted the last animal for this analysis because mnrval keeps producing a NaN probability for one predictor value

% define output variables
tmp_saccade_amp_accuracies = NaN(1,numel(animals));
tmp_endpoint_accuracies = NaN(1,numel(animals));

tmp_saccade_amp_cmtrx = NaN(4,4,numel(animals)); % saccade amp confusion matrix
tmp_endpoint_cmtrx = NaN(4,4,numel(animals)); % saccade endpoint confusion matrix

for AnimalIdx = animals

    mtrx = mtrx_array{AnimalIdx};
    
    % Create stim-evoked filterMtrx(es)
    shank_filters = {'Filters.leftPuff2Filter', 'Filters.leftPuffFilter', 'Filters.rightPuffFilter', 'Filters.rightPuff2Filter'}; % LP2=most anterior site, RP2=most posterior site
    for SHANK = 1:numel(shank_filters)
        trial_param.starts = find(mtrx.allPuff==1);
        [Filters,hankelMtrx] = Create_Filters(mtrx,trial_param,filter_param);
        FILTER = Filters.noNanFilter & Filters.smoothStrainGaugeFilter & eval(shank_filters{SHANK}) & Filters.saccadeAfterAlignFilter;
        filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
        shankFilterMtrx(SHANK) = filterMtrx;
    end
    

    % Calculate  predictor variables for each shank
    for SHANK = 1:num_shanks
        shank(SHANK).init_pos = shankFilterMtrx(SHANK).unalignedPupil(:, trial_param.align);
        shank(SHANK).saccade_amp = find_first_saccade_in_response_window(shankFilterMtrx(SHANK).saccade, trial_param.response_window);
        shank(SHANK).end_pos = shank(SHANK).init_pos + shank(SHANK).saccade_amp;
    end

    % concatenate predictor variables and create response vector
    x1 = [
        shank(1).saccade_amp;
        shank(2).saccade_amp;
        shank(3).saccade_amp;
        shank(4).saccade_amp
        ];

    x2 = [
        shank(1).end_pos;
        shank(2).end_pos;
        shank(3).end_pos;
        shank(4).end_pos
        ];

    y = [
        1*ones(numel(shank(1).saccade_amp),1);
        2*ones(numel(shank(2).saccade_amp),1)
        3*ones(numel(shank(3).saccade_amp),1)
        4*ones(numel(shank(4).saccade_amp),1)
        ];
    
    % z-score standardize predictors
    x1 = (x1-mean(x1))./std(x1);
    x2 = (x2-mean(x2))./std(x2);


    % how many times should 
    num_iterations = 100;

    validation.accuracies.saccade_amp = NaN(num_iterations,1);
    validation.accuracies.endpoint = NaN(num_iterations,1);
    
    validation.confusionmat.saccade_amp = NaN(4,4,num_iterations);
    validation.confusionmat.endpoint = NaN(4,4,num_iterations);
    
    for iteration = 1:num_iterations

        % downsaple majority class
        predictors = [x1, x2];
        response = y;
        [response,predictors] = downsample_majority_class(response, predictors);

        % create cross validation partitions
        KFolds = 5;
        cvp = cvpartition(response, 'KFold', KFolds);

        
        % CROSS-VALIDATE SACCADE AMPLITUDE MODEL
        validationPredictions = NaN(size(response));
        validationScores = NaN(size(response));

        predictor = predictors(:,1);
        for fold = 1:KFolds
            trainingPredictors = predictor(cvp.training(fold), :);
            trainingResponse = response(cvp.training(fold), :);

            cv_mdl = train_multinomial_classifier(trainingPredictors, trainingResponse);

            testPredictors = predictor(cvp.test(fold), :);
            [foldPredictions, foldScores] = cv_mdl.predictFcn(testPredictors);

            % Store predictions in the original order
            validationPredictions(cvp.test(fold), :) = foldPredictions;
        end

        % Compute validation accuracy
        correctPredictions = (validationPredictions == response);
        validationAccuracy = sum(correctPredictions)/length(correctPredictions);
        validationConfusionMatrix = confusionmat(response,validationPredictions);
    %     figure; confusionchart(validationConfusionMatrix)
        validation.accuracies.saccade_amp(iteration) = validationAccuracy;
        validation.confusionmat.saccade_amp(:,:,iteration) = validationConfusionMatrix/sum(validationConfusionMatrix(:));


        % CROSS-VALIDATE ENDPOINT MODEL
        validationPredictions = NaN(size(response));
        validationScores = NaN(size(response));

        predictor = predictors(:,2);
        for fold = 1:KFolds
            trainingPredictors = predictor(cvp.training(fold), :);
            trainingResponse = response(cvp.training(fold), :);

            cv_mdl = train_multinomial_classifier(trainingPredictors, trainingResponse);

            testPredictors = predictor(cvp.test(fold), :);
            [foldPredictions, foldScores] = cv_mdl.predictFcn(testPredictors);

            % Store predictions in the original order
            validationPredictions(cvp.test(fold), :) = foldPredictions;
        end

        % Compute validation accuracy
        correctPredictions = (validationPredictions == response);
        validationAccuracy = sum(correctPredictions)/length(correctPredictions);
        validationConfusionMatrix = confusionmat(response,validationPredictions);
    %     figure; confusionchart(validationConfusionMatrix)
        validation.accuracies.endpoint(iteration) = validationAccuracy;
        validation.confusionmat.endpoint(:,:,iteration) = validationConfusionMatrix/sum(validationConfusionMatrix(:));

    end
    
    tmp_saccade_amp_accuracies(AnimalIdx) = mean(validation.accuracies.saccade_amp);
    tmp_endpoint_accuracies(AnimalIdx) = mean(validation.accuracies.endpoint);
    
    tmp_saccade_amp_cmtrx(:,:,AnimalIdx) = mean(validation.confusionmat.saccade_amp, 3);
    tmp_endpoint_cmtrx(:,:,AnimalIdx) = mean(validation.confusionmat.endpoint, 3);
    
    fprintf('Average model accuracy when using saccade amplitude: %.2f (%.2f)\n', mean(validation.accuracies.saccade_amp), std(validation.accuracies.saccade_amp))
    fprintf('Average model accuracy when using saccade endpoint: %.2f (%.2f)\n', mean(validation.accuracies.endpoint), std(validation.accuracies.endpoint))

end


% PLOT VALIDATION ACCURACIES
figure; set(gcf, 'Position',  [200, 300, 200, 200]); hold on
bar(1,mean(mean(tmp_saccade_amp_accuracies, 2)), 'FaceColor', c_pupil_light)
bar(2,mean(mean(tmp_endpoint_accuracies, 2)), 'FaceColor', c_pupil_light)
h = plot(1:2, [tmp_saccade_amp_accuracies; tmp_endpoint_accuracies], 'o-', 'Color', [0 0 0], 'MarkerSize',4);
xticks([1, 2, 3])
xticklabels({'Amplitude', 'Endpoint'})
xtickangle(25)
ylabel('Accuracy')
ylim([0 1])
hline(0.25, 'k-')


% STATISTICS
fprintf('\n\nSummary Statistics\n')
fprintf('Saccade amplitude shank classifier accuracy : %f +- %f (mean +- SD)\n', mean(tmp_saccade_amp_accuracies), std(tmp_saccade_amp_accuracies))
fprintf('Saccade endpoint shank classifier accuracy : %f +- %f (mean +- SD)\n', mean(tmp_endpoint_accuracies), std(tmp_endpoint_accuracies))

[~, p] = ttest2(tmp_saccade_amp_accuracies, tmp_endpoint_accuracies);
fprintf('t-test2: classifier accuracies: saccade amp vs endpoint p: %f\n', p)

[~, p] = ttest(tmp_saccade_amp_accuracies, 0.25);
fprintf('t-test: classifier accuracies: saccade amp vs random chance (0.25) p: %f\n', p)

[~, p] = ttest(tmp_endpoint_accuracies, 0.25);
fprintf('t-test: classifier accuracies: saccade endpoint vs random chance (0.25) p: %f\n', p)

%% PLOT SACCADE PROBABILITIES

figure; set(gcf, 'Position',  [200, 300, 800, 200]); hold on
for SHANK = 1:num_shanks
    

    greens = cbrewer('seq','YlGn',4, 'PCHIP');
    magentas = cbrewer('seq','RdPu',4, 'PCHIP');
    left_color = greens(3,:);
    right_color = magentas(3,:);

    % set NaNs to zero
    shank(SHANK).probabilities.saccade_left(isnan(shank(SHANK).probabilities.saccade_left)) = 0;
    shank(SHANK).probabilities.saccade_right(isnan(shank(SHANK).probabilities.saccade_right)) = 0;
    
    % Plot saccade probability
    subplot(1,4,SHANK); hold on
    plot(shank(SHANK).probabilities.saccade_bin_centers', shank(SHANK).probabilities.saccade_left', 'Color', left_color)
    plot(mean(shank(SHANK).probabilities.saccade_bin_centers', 2), mean(shank(SHANK).probabilities.saccade_left', 2), 'Color', left_color, 'LineWidth', 3)
    
    plot(shank(SHANK).probabilities.saccade_bin_centers', shank(SHANK).probabilities.saccade_right', 'Color', right_color)
    plot(mean(shank(SHANK).probabilities.saccade_bin_centers', 2), mean(shank(SHANK).probabilities.saccade_right', 2), 'Color', right_color, 'LineWidth', 3)
    ylim([0 1])
    ylabel('P(Saccade)')
    xlabel('Initial pupil angle (°)')
    axis square

end

%% PLOT HEAD MOVEMENT PROBABILITIES

figure; set(gcf, 'Position',  [200, 300, 800, 200]); hold on
for SHANK = 1:num_shanks
    

    greens = cbrewer('seq','YlGn',4, 'PCHIP');
    magentas = cbrewer('seq','RdPu',4, 'PCHIP');
    left_color = greens(3,:);
    right_color = magentas(3,:);

    % set NaNs to zero
    shank(SHANK).probabilities.head_mov_left(isnan(shank(SHANK).probabilities.head_mov_left)) = 0;
    shank(SHANK).probabilities.head_mov_right(isnan(shank(SHANK).probabilities.head_mov_right)) = 0;
    
    % Plot saccade probability
    subplot(1,4,SHANK); hold on
    plot(shank(SHANK).probabilities.head_mov_bin_centers', shank(SHANK).probabilities.head_mov_left', 'Color', left_color)
    plot(mean(shank(SHANK).probabilities.head_mov_bin_centers', 2), mean(shank(SHANK).probabilities.head_mov_left', 2), 'Color', left_color, 'LineWidth', 3)
    
    plot(shank(SHANK).probabilities.head_mov_bin_centers', shank(SHANK).probabilities.head_mov_right', 'Color', right_color)
    plot(mean(shank(SHANK).probabilities.head_mov_bin_centers', 2), mean(shank(SHANK).probabilities.head_mov_right', 2), 'Color', right_color, 'LineWidth', 3)
    ylim([0 1])
    ylabel('P(Head movement)')
    xlabel('Initial pupil angle (°)')
    axis square

end
