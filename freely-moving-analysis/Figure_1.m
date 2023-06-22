%% Load data
clear all

BehaviorFilePath = '../Data/Figure 1 and supplements/Freely moving SC stimulation';
Experiment = {

    'SZ399', [20220425:20220427], [01]; % SC ChR2
    'SZ409', [20220504:20220506], [01]; % SC ChR2
    'SZ410', [20220504:20220506], [01]; % SC ChR2
    'SZ413', [20220504:20220506], [01]; % SC ChR2
    'SZ415', [20220504:20220506], [01]; % SC ChR2 (EXAMPLE ANIMAL)
    'SZ416', [20220504:20220506], [01]; % SC ChR2
    'SZ417', [20220504:20220506], [01]; % SC ChR2

    };

% Set trial parameters for event detection (these are used by combine_sessions_freely_moving)
trial_param.saccade_threshold = 2; % deg
trial_param.velocity_threshold = 550; % deg/sec
trial_param.saccade_onset_velocity_threshold = 70; % deg/sec
trial_param.saccade_offset_velocity_threshold = 60; % deg/sec
trial_param.head_movement_velocity_threshold = 140; % deg/sec
trial_param.head_movement_minimum_duration = 3; % frames
trial_param.frame_rate = 100; % Hz (default: 100)

% Set trial parameters for trial structure
trial_param.window_size = 61; % Num frames (width of each trial for analysis)
trial_param.align = 11; % Frame between 1 and window_size. Specifies which frame you want to define as 0
trial_param.response_window = trial_param.align:trial_param.align+10; % Frame range for "response window," which is important for filters
trial_param.preAlignPeriod = 1:trial_param.align-1;
trial_param.window_xaxis = ([1:trial_param.window_size] - trial_param.align)*1/trial_param.frame_rate;

% Load data and create finalized time series
for AnimalIdx = 1:size(Experiment, 1)
    DateRange = [Experiment{AnimalIdx,2}];
    AnimalList = Experiment{AnimalIdx,1};
    SessionNumber = Experiment{AnimalIdx,3};
    BehaviorFiles=LoadBehavior('BehaviorFilePath', BehaviorFilePath, 'AnimalList', AnimalList, 'DateRange', DateRange,'SessionNumber',SessionNumber);
    mtrx_array{AnimalIdx} = combine_sessions_freely_moving(BehaviorFiles, trial_param);
end

%% PLOT EXAMPLE TRACES
% Example used: Animal SZ399, filtered trial 12 (trial 154 unfiltered)
AnimalIdx = 1;
TRIAL = 12;

% filter parameters
filter_param = struct();
filter_param.head_movement_threshold = 1;
filter_param.smooth_head_window = trial_param.align-10:trial_param.align;
filter_param.pupil_movement_threshold = 6;
filter_param.smooth_pupil_window = trial_param.align-5:trial_param.align;
    
% Select mtrx for single animal
mtrx = mtrx_array{AnimalIdx};

% Create filterMtrx (opto trials aligned to trial start)
trial_param.starts = find(mtrx.optoTrial==1);
[Filters,hankelMtrx] = Create_Filters(mtrx,trial_param,filter_param);
FILTER = Filters.noNanFilter & Filters.noMovementBeforeFirstSaccadeFilter & Filters.simultaneousSaccadeAfterAlignFilter & Filters.smoothPupilFilter & Filters.smoothHeadFilter & Filters.headYawEventAfterAlignFilter;
filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);

% FORMATTING
head_color = [1, 0.7, 0];
eye_color = [0.5 0 0.5];
gaze_color = [0.5 0.5 0.5];

figure; set(gcf, 'Position',  [200, 100, 300, 300]); hold on
head_trace = filterMtrx.headBodyYawAngleIMU(TRIAL,:)- filterMtrx.headBodyYawAngleIMU(TRIAL,trial_param.align) + filterMtrx.headBodyYawAngleVid(TRIAL,trial_param.align);
pupil_trace = filterMtrx.meanPupilX(TRIAL,:);
gaze_trace = head_trace + pupil_trace;

plot(trial_param.window_xaxis, head_trace, 'Color', head_color);
xlim([trial_param.window_xaxis(1), trial_param.window_xaxis(end)])
ylim([-60, 60])

plot(trial_param.window_xaxis, pupil_trace, 'Color', eye_color)
xlim([trial_param.window_xaxis(1), trial_param.window_xaxis(end)])
ylim([-60, 60])

plot(trial_param.window_xaxis, gaze_trace, 'Color', gaze_color)
xlim([trial_param.window_xaxis(1), trial_param.window_xaxis(end)])
ylim([-60, 60])

vline(0)
vline(0.1)


%% SACCADE AND HEAD PROBABILITIES

% Set filter paramters
filter_param = struct();
filter_param.head_movement_threshold = 1;
filter_param.smooth_head_window = trial_param.align-10:trial_param.align;
filter_param.pupil_movement_threshold = 6;
filter_param.smooth_pupil_window = trial_param.align-5:trial_param.align;

% Formatting
greens = cbrewer('seq','YlGn',4, 'PCHIP');
magentas = cbrewer('seq','RdPu',4, 'PCHIP');
left_color = greens(3,:);
right_color = magentas(3,:);

% Define bin edges
num_bins = 4;
pupil_bin_range = [-60, 60];
pupil_bin_edges = linspace(pupil_bin_range(1), pupil_bin_range(2), num_bins+1);
head_bin_edges = pupil_bin_edges;

% Initialize aggregate variables
probabilities = [];

for AnimalIdx = 1:7
    mtrx = mtrx_array{AnimalIdx};
    
    % Create filterMtrx (aligned to trial start)
    trial_param.starts = find(mtrx.optoTrial==1);
    [Filters,hankelMtrx] = Create_Filters(mtrx,trial_param,filter_param);
    
    % SACCADE PROBABILITIES: get counts for all trials and bin centers
    FILTER = Filters.noNanFilter & Filters.smoothPupilFilter & Filters.smoothHeadFilter;
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    [N_pupil_all, pupil_bin_edges] = histcounts(filterMtrx.meanPupilX(:,trial_param.align), pupil_bin_edges);
    pupil_bin_width = diff(pupil_bin_edges);
    pupil_bin_width = pupil_bin_width(1);
    pupil_bin_centers = pupil_bin_edges(1:end-1)+pupil_bin_width/2;
    probabilities.saccade_bin_centers(AnimalIdx, :) = pupil_bin_centers;
    
    % all saccades
    FILTER = Filters.noNanFilter & Filters.noMovementBeforeFirstSaccadeFilter & Filters.simultaneousSaccadeAfterAlignFilter & Filters.smoothPupilFilter & Filters.smoothHeadFilter;
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    N_saccade = histcounts(filterMtrx.meanPupilX(:,trial_param.align), pupil_bin_edges);
    probabilities.saccade(AnimalIdx, :) = N_saccade./N_pupil_all;

    % left saccades
    FILTER = Filters.noNanFilter & Filters.noMovementBeforeFirstSaccadeFilter & Filters.simultaneousSaccadeAfterAlignFilter & Filters.bothPupilLeftSaccadeAfterAlignFilter & Filters.smoothPupilFilter & Filters.smoothHeadFilter;
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    N_saccade = histcounts(filterMtrx.meanPupilX(:,trial_param.align), pupil_bin_edges);
    probabilities.saccade_left(AnimalIdx, :) = N_saccade./N_pupil_all;

    % right saccades
    FILTER = Filters.noNanFilter & Filters.noMovementBeforeFirstSaccadeFilter & Filters.simultaneousSaccadeAfterAlignFilter & Filters.bothPupilRightSaccadeAfterAlignFilter & Filters.smoothPupilFilter & Filters.smoothHeadFilter;
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    N_saccade = histcounts(filterMtrx.meanPupilX(:,trial_param.align), pupil_bin_edges);
    probabilities.saccade_right(AnimalIdx, :) = N_saccade./N_pupil_all;
    

    % HEAD PROBABILITIES: get counts for all trials and bin centers
    FILTER = Filters.noNanFilter & Filters.smoothPupilFilter & Filters.smoothHeadFilter;
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    [N_head_all, head_bin_edges] = histcounts(filterMtrx.headBodyYawAngleIMU(:,trial_param.align), head_bin_edges);
    head_bin_width = diff(head_bin_edges);
    head_bin_width = head_bin_width(1);
    head_bin_centers = head_bin_edges(1:end-1)+head_bin_width/2;
    probabilities.head_bin_centers(AnimalIdx, :) = head_bin_centers;
    
    % all head movements
    FILTER = Filters.noNanFilter & Filters.smoothPupilFilter & Filters.smoothHeadFilter & Filters.headYawEventAfterAlignFilter;
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    N_head = histcounts(filterMtrx.headBodyYawAngleIMU(:,trial_param.align), head_bin_edges);
    probabilities.head(AnimalIdx, :) = N_head./N_head_all;

    % left head movement
    FILTER = Filters.noNanFilter & Filters.smoothPupilFilter & Filters.smoothHeadFilter & Filters.leftHeadYawEventAfterAlignFilter;
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    N_head = histcounts(filterMtrx.headBodyYawAngleIMU(:,trial_param.align), head_bin_edges);
    probabilities.head_left(AnimalIdx, :) = N_head./N_head_all;

    % right head movement
    FILTER = Filters.noNanFilter & Filters.smoothPupilFilter & Filters.smoothHeadFilter & Filters.rightHeadYawEventAfterAlignFilter;
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);
    N_head = histcounts(filterMtrx.headBodyYawAngleIMU(:,trial_param.align), head_bin_edges);
    probabilities.head_right(AnimalIdx, :) = N_head./N_head_all;
    
end

% set NaNs to zero
probabilities.saccade_left(isnan(probabilities.saccade_left)) = 0;
probabilities.saccade_right(isnan(probabilities.saccade_right)) = 0;
probabilities.head_left(isnan(probabilities.head_left)) = 0;
probabilities.head_right(isnan(probabilities.head_right)) = 0;

% Plot saccade probability
figure; set(gcf, 'Position',  [200, 100, 250, 250]); hold on
plot(probabilities.saccade_bin_centers', probabilities.saccade_left', 'Color', left_color)
plot(mean(probabilities.saccade_bin_centers', 2), mean(probabilities.saccade_left', 2), 'Color', left_color, 'LineWidth', 3)

plot(probabilities.saccade_bin_centers', probabilities.saccade_right', 'Color', right_color)
plot(mean(probabilities.saccade_bin_centers', 2), mean(probabilities.saccade_right', 2), 'Color', right_color, 'LineWidth', 3)
ylim([0 1])
ylabel('P(Saccade)')
xlabel('Initial pupil angle (°)')

% Plot head movement probability
figure; set(gcf, 'Position',  [200, 100, 250, 250]); hold on
plot(probabilities.head_bin_centers', probabilities.head_left', 'Color', left_color)
plot(mean(probabilities.head_bin_centers', 2), mean(probabilities.head_left', 2), 'Color', left_color, 'LineWidth', 3)

plot(probabilities.head_bin_centers', probabilities.head_right', 'Color', right_color)
plot(mean(probabilities.head_bin_centers', 2), mean(probabilities.head_right', 2), 'Color', right_color, 'LineWidth', 3)
ylim([0 1])
ylabel('P(Head mov.)')
xlabel('Initial head angle (°)')

%% GAZE/PUPIL/HEAD MOVEMENT ANALYSIS

% Specify which mice to make individual plots
% plot_mice = 1:numel(mtrx_array); % plot all mice
plot_mice = 5; % example mouse

% filter parameters
filter_param = struct();
filter_param.head_movement_threshold = 1;
filter_param.smooth_head_window = trial_param.align-10:trial_param.align;
filter_param.pupil_movement_threshold = 6;
filter_param.smooth_pupil_window = trial_param.align-5:trial_param.align;

% initialize aggregate variables
regression_stats = [];
population_averages = [];

for AnimalIdx = 1:numel(mtrx_array)
    
    % Select mtrx for single animal
    mtrx = mtrx_array{AnimalIdx};

    % Create filterMtrx (opto trials aligned to trial start)
    trial_param.starts = find(mtrx.optoTrial==1);
    [Filters,hankelMtrx] = Create_Filters(mtrx,trial_param,filter_param);
    FILTER = Filters.noNanFilter & Filters.noMovementBeforeFirstSaccadeFilter & Filters.simultaneousSaccadeAfterAlignFilter & Filters.smoothPupilFilter & Filters.smoothHeadFilter & Filters.headYawEventAfterAlignFilter;
    filterMtrx = Create_filterMtrx(hankelMtrx,FILTER);

    % FEATURE ENGINEERING (i.e., identify movement amplitudes, start positions, etc.)
    
    headStartAngle = filterMtrx.headBodyYawAngleIMU(:,trial_param.align);
    headAmplitude = row_indexing(filterMtrx.headYawEventsAmp, first_event_after(filterMtrx.headYawEventsAmp, trial_param.align));
    headEndAngle = headStartAngle + headAmplitude;
    headEventStartingFrame = first_event_after(filterMtrx.headYawEventsDuration, trial_param.align); % find frame where head movement begins
    headEventDuration = row_indexing(filterMtrx.headYawEventsDuration, first_event_after(filterMtrx.headYawEventsDuration, trial_param.align)); % find duration of head movement (in frames)
    headEndFrame = headEventStartingFrame + headEventDuration;
    
    % head pitch variables
    headPitchStartAngle = filterMtrx.headBodyPitchAngleGyro(:,trial_param.align);
    headPitchEndAngle = row_indexing(filterMtrx.headBodyPitchAngleGyro, headEndFrame);
    headPitchAmplitude = headPitchEndAngle - headPitchStartAngle;
    
    % head roll variables
    headRollStartAngle = filterMtrx.headBodyRollAngleGyro(:,trial_param.align);
    headRollEndAngle = row_indexing(filterMtrx.headBodyRollAngleGyro, headEndFrame);
    headRollAmplitude = headRollEndAngle - headRollStartAngle;

    % left pupil variables
    leftPupilStartAngle = filterMtrx.leftPupilX(:,trial_param.align);
    leftPupilSaccadeAmplitude = row_indexing(filterMtrx.leftPupilXSaccades, first_event_after(filterMtrx.leftPupilXSaccades, trial_param.align));
    leftPupilSaccadeEndAngle = leftPupilStartAngle + leftPupilSaccadeAmplitude;
    leftPupilVorEndAngle = row_indexing(filterMtrx.leftPupilX, headEndFrame);
    leftPupilSaccadeStartFrame = first_event_after(filterMtrx.leftPupilXSaccades, trial_param.align); % find frame where head movement begins
    leftPupilSaccadeDuration = row_indexing(filterMtrx.leftPupilSaccadesLength, first_event_after(filterMtrx.leftPupilSaccadesLength, trial_param.align)); % find duration of head movement (in frames)
    leftPupilSaccadeEndFrame = leftPupilSaccadeStartFrame + leftPupilSaccadeDuration;
    leftPupilStartAngleY = filterMtrx.leftPupilY(:,trial_param.align);
    leftPupilSaccadeAmplitudeY = row_indexing(filterMtrx.leftPupilYSaccades, first_event_after(filterMtrx.leftPupilYSaccades, trial_param.align));
    
    % right pupil variables
    rightPupilStartAngle = filterMtrx.rightPupilX(:,trial_param.align);
    rightPupilSaccadeAmplitude = row_indexing(filterMtrx.rightPupilXSaccades, first_event_after(filterMtrx.rightPupilXSaccades, trial_param.align));
    rightPupilSaccadeEndAngle = rightPupilStartAngle + rightPupilSaccadeAmplitude;
    rightPupilVorEndAngle = row_indexing(filterMtrx.rightPupilX, headEndFrame);
    rightPupilSaccadeStartFrame = first_event_after(filterMtrx.rightPupilXSaccades, trial_param.align); % find frame where head movement begins
    rightPupilSaccadeDuration = row_indexing(filterMtrx.rightPupilSaccadesLength, first_event_after(filterMtrx.rightPupilSaccadesLength, trial_param.align)); % find duration of head movement (in frames)
    rightPupilSaccadeEndFrame = rightPupilSaccadeStartFrame + rightPupilSaccadeDuration;
    rightPupilStartAngleY = filterMtrx.rightPupilX(:,trial_param.align);
    rightPupilSaccadeAmplitudeY = row_indexing(filterMtrx.rightPupilYSaccades, first_event_after(filterMtrx.rightPupilYSaccades, trial_param.align));
    
    % mean pupil variables
    meanPupilSaccadeStartFrame = round((leftPupilSaccadeStartFrame + rightPupilSaccadeStartFrame)/2);
    meanPupilSaccadeDuration = round((leftPupilSaccadeDuration + rightPupilSaccadeDuration)/2);
    meanPupilSaccadeEndFrame = round((leftPupilSaccadeEndFrame + rightPupilSaccadeEndFrame)/2);
    meanPupilStartAngle = row_indexing(filterMtrx.meanPupilX, meanPupilSaccadeStartFrame);
    meanPupilSaccadeEndAngle = row_indexing(filterMtrx.meanPupilX, meanPupilSaccadeEndFrame);
    meanPupilSaccadeAmplitude = meanPupilSaccadeEndAngle - meanPupilStartAngle;
    meanPupilVorEndAngle = row_indexing(filterMtrx.meanPupilX, headEndFrame);
    meanPupilOverallAmplitude = meanPupilVorEndAngle - meanPupilStartAngle;
    meanPupilStartAngleY = row_indexing(filterMtrx.meanPupilY, meanPupilSaccadeStartFrame);
    meanPupilSaccadeEndAngleY = row_indexing(filterMtrx.meanPupilY, meanPupilSaccadeEndFrame);
    meanPupilSaccadeAmplitudeY = meanPupilSaccadeEndAngleY - meanPupilStartAngleY;
    
    % mean pupil velocity variables
    meanPupilXVelocity = [];
    meanPupilMaxSacadeVelocity = [];
    for i = 1:size(filterMtrx.meanPupilX,1)
        meanPupilXVelocity(i,:) = [0, diff(filterMtrx.meanPupilX(i,:))]/(1/trial_param.frame_rate); 
        meanPupilMaxSacadeVelocity(i,:) = max(abs(meanPupilXVelocity(i, meanPupilSaccadeStartFrame(i):meanPupilSaccadeEndFrame(i))));
    end

    gazeHeadComponentAmplitude = row_indexing(filterMtrx.headBodyYawAngleIMU, meanPupilSaccadeEndFrame) - headStartAngle;
    gazeStartAngle = headStartAngle + meanPupilStartAngle;
    gazeAmplitude = gazeHeadComponentAmplitude + meanPupilSaccadeAmplitude;
    gazeSaccadeEndAngle = row_indexing(filterMtrx.headBodyYawAngleIMU, meanPupilSaccadeEndFrame) + meanPupilSaccadeEndAngle;
    gazeVorEndAngle = headEndAngle + meanPupilVorEndAngle; % endpoint of head+meanPupil at the frame when the head movement stops

    % FORMATTING
    c_pupil = [127, 63, 152]/255;
    c_head = [247, 148, 29]/255;
    c_gaze = [128, 128, 128]/255;
    c_pupil_light = (1 - c_pupil)*0.5 + c_pupil;
    c_head_light = (1 - c_head)*0.5 + c_head;
    c_gaze_light = (1 - c_gaze)*0.5 + c_gaze;

    % MOVEMENT TRACE PLOTS
    if any(AnimalIdx == plot_mice)
        x_lims = [-0.1 0.2];
        y_lims = [-80 80];

        % gaze traces
        figure; hold on; set(gcf, 'Position',  [200, 100, 200, 200]); hold on
        plot(trial_param.window_xaxis, (filterMtrx.meanPupilX + filterMtrx.headBodyYawAngleIMU)', 'Color', c_gaze)
        xlim(x_lims); ylim(y_lims)
        vline(0, 'k--')
        xlabel('Time (s)')
        ylabel('Gaze angle (°)')

        % head traces
        figure; hold on; set(gcf, 'Position',  [200, 100, 200, 200]); hold on
        plot(trial_param.window_xaxis, filterMtrx.headBodyYawAngleIMU', 'Color', c_head)
        xlim(x_lims); ylim(y_lims)
        vline(0, 'k--')
        xlabel('Time (s)')
        ylabel('Head-body angle (°)')

        % pupil traces
        figure; hold on; set(gcf, 'Position',  [200, 100, 200, 200]); hold on
        plot(trial_param.window_xaxis, filterMtrx.meanPupilX', 'Color', c_pupil)
        xlim(x_lims); ylim(y_lims)
        vline(0, 'k--')
        xlabel('Time (s)')
        ylabel('Pupil angle (°)')
    end

    % SPECIFY PUPIL VARIABLES
    % Change to plot mean, left, or right pupil (make sure all are the same)
    PupilStartAngleForPlotting = meanPupilStartAngle; 
    PupilSaccadeAmplitudeForPlotting = meanPupilSaccadeAmplitude;

    % GENERATE SCATTER PLOTS FOR INDIVIDUAL MICE
    if any(AnimalIdx == plot_mice)

        % Limits for most scatter plots
        xlimits = [-100, 100];
        ylimits = [-100, 100];

        % MAIN FIGURES
        scatter_plot_with_regression(gazeStartAngle, gazeAmplitude, 'Color', c_gaze, 'XLim', xlimits, 'YLim', ylimits);
        xlabel('Initial gaze angle (°)'); ylabel('Gaze displacement (°)')

        scatter_plot_with_regression(headStartAngle, gazeAmplitude, 'Color', c_gaze, 'XLim', xlimits, 'YLim', ylimits);
        xlabel('Initial head angle (°)'); ylabel('Gaze displacement (°)')

        scatter_plot_with_regression(PupilStartAngleForPlotting, gazeAmplitude, 'Color', c_gaze, 'XLim', xlimits, 'YLim', ylimits);
        xlabel('Initial pupil angle (°)'); ylabel('Gaze displacement (°)')


        % NOTE: using gazeHeadComponentAmplitude instead of headAmplitude
        scatter_plot_with_regression(gazeStartAngle, gazeHeadComponentAmplitude, 'Color', c_head, 'XLim', xlimits, 'YLim', ylimits);
        xlabel('Initial gaze angle (°)'); ylabel('Head displacement (°)')

        scatter_plot_with_regression(headStartAngle, gazeHeadComponentAmplitude, 'Color', c_head, 'XLim', xlimits, 'YLim', ylimits);
        xlabel('Initial head angle (°)'); ylabel('Head displacement (°)')

        scatter_plot_with_regression(PupilStartAngleForPlotting, gazeHeadComponentAmplitude, 'Color', c_head, 'XLim', xlimits, 'YLim', ylimits);
        xlabel('Initial pupil angle (°)'); ylabel('Head displacement (°)')

        
        scatter_plot_with_regression(gazeStartAngle, PupilSaccadeAmplitudeForPlotting, 'Color', c_pupil, 'XLim', xlimits, 'YLim', ylimits);
        xlabel('Initial gaze angle (°)'); ylabel('Saccade displacement (°)')

        scatter_plot_with_regression(headStartAngle, PupilSaccadeAmplitudeForPlotting, 'Color', c_pupil, 'XLim', xlimits, 'YLim', ylimits);
        xlabel('Initial head angle (°)'); ylabel('Saccade displacement (°)')

        scatter_plot_with_regression(PupilStartAngleForPlotting, PupilSaccadeAmplitudeForPlotting, 'Color', c_pupil, 'XLim', xlimits, 'YLim', ylimits);
        xlabel('Initial pupil angle (°)'); ylabel('Saccade displacement (°)')

    end

    % CALCULATE REGRESSION R2 VALUES

    mdl = fitlm(gazeStartAngle, gazeAmplitude);
    regression_stats.gaze_v_gaze.r2(AnimalIdx,1) = mdl.Rsquared.Ordinary;
    regression_stats.gaze_v_gaze.coef(AnimalIdx,1) = mdl.Coefficients.Estimate(2);
    regression_stats.gaze_v_gaze.p(AnimalIdx,1) = mdl.Coefficients.pValue(2);
    
    mdl = fitlm(PupilStartAngleForPlotting, gazeAmplitude);
    regression_stats.pupil_v_gaze.r2(AnimalIdx,1) = mdl.Rsquared.Ordinary;
    regression_stats.pupil_v_gaze.coef(AnimalIdx,1) = mdl.Coefficients.Estimate(2);
    regression_stats.pupil_v_gaze.p(AnimalIdx,1) = mdl.Coefficients.pValue(2);

    mdl = fitlm(headStartAngle, gazeAmplitude);
    regression_stats.head_v_gaze.r2(AnimalIdx,1) = mdl.Rsquared.Ordinary;
    regression_stats.head_v_gaze.coef(AnimalIdx,1) = mdl.Coefficients.Estimate(2);
    regression_stats.head_v_gaze.p(AnimalIdx,1) = mdl.Coefficients.pValue(2);

    % NOTE: using gazeHeadComponentAmplitude instead of headAmplitude
    mdl = fitlm(gazeStartAngle, gazeHeadComponentAmplitude);
    regression_stats.gaze_v_head.r2(AnimalIdx,1) = mdl.Rsquared.Ordinary;
    regression_stats.gaze_v_head.coef(AnimalIdx,1) = mdl.Coefficients.Estimate(2);
    regression_stats.gaze_v_head.p(AnimalIdx,1) = mdl.Coefficients.pValue(2);
    
    mdl = fitlm(PupilStartAngleForPlotting, gazeHeadComponentAmplitude);
    regression_stats.pupil_v_head.r2(AnimalIdx,1) = mdl.Rsquared.Ordinary;
    regression_stats.pupil_v_head.coef(AnimalIdx,1) = mdl.Coefficients.Estimate(2);
    regression_stats.pupil_v_head.p(AnimalIdx,1) = mdl.Coefficients.pValue(2);

    mdl = fitlm(headStartAngle, gazeHeadComponentAmplitude);
    regression_stats.head_v_head.r2(AnimalIdx,1) = mdl.Rsquared.Ordinary;
    regression_stats.head_v_head.coef(AnimalIdx,1) = mdl.Coefficients.Estimate(2);
    regression_stats.head_v_head.p(AnimalIdx,1) = mdl.Coefficients.pValue(2);


    mdl = fitlm(gazeStartAngle, PupilSaccadeAmplitudeForPlotting);
    regression_stats.gaze_v_pupil.r2(AnimalIdx,1) = mdl.Rsquared.Ordinary;
    regression_stats.gaze_v_pupil.coef(AnimalIdx,1) = mdl.Coefficients.Estimate(2);
    regression_stats.gaze_v_pupil.p(AnimalIdx,1) = mdl.Coefficients.pValue(2);

    mdl = fitlm(PupilStartAngleForPlotting, PupilSaccadeAmplitudeForPlotting);
    regression_stats.pupil_v_pupil.r2(AnimalIdx,1) = mdl.Rsquared.Ordinary;
    regression_stats.pupil_v_pupil.coef(AnimalIdx,1) = mdl.Coefficients.Estimate(2);
    regression_stats.pupil_v_pupil.p(AnimalIdx,1) = mdl.Coefficients.pValue(2);

    mdl = fitlm(headStartAngle, PupilSaccadeAmplitudeForPlotting);
    regression_stats.head_v_pupil.r2(AnimalIdx,1) = mdl.Rsquared.Ordinary;
    regression_stats.head_v_pupil.coef(AnimalIdx,1) = mdl.Coefficients.Estimate(2);
    regression_stats.head_v_pupil.p(AnimalIdx,1) = mdl.Coefficients.pValue(2);  


    % SUPPLEMENT: SACCADE HORIZONTAL/VERTICAL SACCADE MAGNITUDES
    % Collect aggregates (figure plotted later)
    left_pupil_saccade_amp_mag_x_mean(AnimalIdx,:) = mean(abs(leftPupilSaccadeAmplitude));
    left_pupil_saccade_amp_mag_y_mean(AnimalIdx,:) = mean(abs(leftPupilSaccadeAmplitudeY));

    right_pupil_saccade_amp_mag_x_mean(AnimalIdx,:) = mean(abs(rightPupilSaccadeAmplitude));
    right_pupil_saccade_amp_mag_y_mean(AnimalIdx,:) = mean(abs(rightPupilSaccadeAmplitudeY));


    % SUPPLEMENT: INTEROCULAR SACCADE CORRELATION
    % Collect aggregates (figure plotted later)
    R = corrcoef(leftPupilSaccadeAmplitude, rightPupilSaccadeAmplitude);
    pupil_saccade_x_correlation(AnimalIdx,:) = R(2,1);
    R = corrcoef(leftPupilSaccadeAmplitudeY, rightPupilSaccadeAmplitudeY);
    pupil_saccade_y_correlation(AnimalIdx,:) = R(2,1);
    

    % SUPPLEMENT: HEAD YAW/ROLL/PITCH MAGNITUDE
    % Collect aggregates (figure plotted later)
    head_yaw_amp_mag_mean(AnimalIdx,:) = mean(abs(headAmplitude));
    head_roll_amp_mag_mean(AnimalIdx,:) = mean(abs(headRollAmplitude));
    head_pitch_amp_mag_mean(AnimalIdx,:) = mean(abs(headPitchAmplitude));


    % SUPPLEMENT: MAIN SEQUENCE

    % Calculate regression values (figure plotted later)
    mdl = fitlm(abs(meanPupilSaccadeAmplitude), abs(meanPupilMaxSacadeVelocity));
    saccade_amp_vs_peak_velocity.r2(AnimalIdx,:) = mdl.Rsquared.Ordinary;
    saccade_amp_vs_peak_velocity.coef(AnimalIdx,:) = mdl.Coefficients.Estimate(2);

    % Scatter plot for individual mice
    % (Note: SZ417 was used as the example main sequence)
    if any(AnimalIdx == plot_mice)
        scatter_plot_with_regression(abs(meanPupilSaccadeAmplitude), abs(meanPupilMaxSacadeVelocity), 'Color', [0 0 0], 'XLim', [0, 60], 'YLim', [0, 3000], 'PlotAxes', '');
        xlabel('Saccade displacement (°)'); ylabel('Saccade velocity (°)')
    end


    % SUPPLEMENT: INITIAL HEAD VS INITIAL PUPIL
    
    % Calculate regression values (figure plotted later)
    mdl = fitlm(headStartAngle, meanPupilStartAngle);
    initial_head_vs_initial_pupil.r2(AnimalIdx,:) = mdl.Rsquared.Ordinary;
    initial_head_vs_initial_pupil.coef(AnimalIdx,:) = mdl.Coefficients.Estimate(2);

    % Scatter plot for individual mice
    if any(AnimalIdx == plot_mice)
        scatter_plot_with_regression(headStartAngle, PupilStartAngleForPlotting, 'Color', [0 0 0], 'XLim', xlimits, 'YLim', ylimits);
        xlabel('Initial head angle (°)'); ylabel('Initial pupil angle (°)')
    end


    % SUPPLEMENT: INITIAL AND FINAL POSITIONS FOR POPULATION (GAZE, HEAD, PUPIL)

    population_averages.initial_gaze(AnimalIdx,:) = mean(gazeStartAngle);
    population_averages.initial_head(AnimalIdx,:) = mean(headStartAngle);
    population_averages.initial_eye(AnimalIdx,:) = mean(PupilStartAngleForPlotting);

    population_averages.final_gaze(AnimalIdx,:) = mean(gazeStartAngle + gazeAmplitude);
    population_averages.final_head(AnimalIdx,:) = mean(headStartAngle + gazeHeadComponentAmplitude);
    population_averages.final_eye(AnimalIdx,:) = mean(PupilStartAngleForPlotting + PupilSaccadeAmplitudeForPlotting);

    population_averages.amplitude_gaze(AnimalIdx,:) = mean(gazeAmplitude);
    population_averages.amplitude_head(AnimalIdx,:) = mean(gazeHeadComponentAmplitude);
    population_averages.amplitude_eye(AnimalIdx,:) = mean(PupilSaccadeAmplitudeForPlotting);


    % SUPPLEMENT: SCATTERS OF INITIAL POSITIONS AND AMPLITUDES
    xlimits = [-100, 100];
    ylimits = [-100, 100];

    figure; set(gcf, 'Position',  [500, 700, 200, 400]); hold on
    
    ax1 = subplot(2,1,1); hold on
    scatter_plot_with_regression(headStartAngle, gazeHeadComponentAmplitude, 'Color', c_head, 'XLim', xlimits, 'YLim', ylimits, 'ax', ax1);
    ylabel('Head displacement (°)')
    xlabel('Initial head angle (°)')

    ax2 = subplot(2,1,2); hold on
    scatter_plot_with_regression(PupilStartAngleForPlotting, PupilSaccadeAmplitudeForPlotting, 'Color', c_pupil, 'XLim', xlimits, 'YLim', ylimits, 'ax', ax2);
    ylabel('Saccade displacement (°)')
    xlabel('Initial pupil angle (°)')

    % SUPPLEMENT: HISTOGRAMS OF INITIAL AND FINAL PUPIL POSITIONS
    figure; set(gcf, 'Position',  [500, 500, 200, 200]); hold on
    bin_edges = [-60:15:60];
    
    subplot(2,1,1); hold on
    histogram(headStartAngle, bin_edges, 'DisplayStyle', 'stairs', 'EdgeColor', c_head, 'LineWidth', 1.5)
    histogram(headStartAngle + gazeHeadComponentAmplitude, bin_edges, 'FaceColor', c_head_light, 'LineStyle', 'none')
    ylabel('count')
    xlabel('Position (°)')

    subplot(2,1,2); hold on
    histogram(PupilStartAngleForPlotting, bin_edges, 'DisplayStyle', 'stairs', 'EdgeColor', c_pupil, 'LineWidth', 1.5)
    histogram(PupilStartAngleForPlotting + PupilSaccadeAmplitudeForPlotting, bin_edges, 'FaceColor', c_pupil_light, 'LineStyle', 'none')
    ylabel('count')
    xlabel('Position (°)')

end

%% PLOT SUMMARY R2 AND SLOPES FOR INITIAL GAZE/HEAD/EYE DEPENDENCE

% Gaze amplitude regressed on initial gaze/head/pupil
tmp_r2_data = [regression_stats.gaze_v_gaze.r2, regression_stats.head_v_gaze.r2, regression_stats.pupil_v_gaze.r2];
tmp_coef_data = [regression_stats.gaze_v_gaze.coef, regression_stats.head_v_gaze.coef, regression_stats.pupil_v_gaze.coef];

figure; set(gcf, 'Position',  [200, 100, 200, 300]); hold on
bar(1:3, mean(tmp_r2_data,1), 'FaceColor', c_gaze_light)
scatter(1:3, tmp_r2_data, 'k')
plot(1:3, tmp_r2_data, 'k-')
ylim([-0.1 1])
xlim([0.5 3.5])
title('Gaze regressed on')
ylabel('R2')
xticks([1 2 3])
xticklabels({'Initial gaze', 'Initial head', 'Intial pupil'})
xtickangle(45)

figure; set(gcf, 'Position',  [200, 100, 200, 300]); hold on
bar(1:3, mean(tmp_coef_data,1), 'FaceColor', c_gaze_light)
scatter(1:3, tmp_coef_data, 'k')
plot(1:3, tmp_coef_data, 'Color', [0.2, 0.2, 0.2])
ylim([-1 1])
xlim([0.5 3.5])
title('Gaze regressed on')
ylabel('Slope')
xticks([1 2 3])
xticklabels({'Initial gaze', 'Initial head', 'Intial pupil'})
xtickangle(45)


% Head amplitude regressed on initial gaze/head/pupil
tmp_r2_data = [regression_stats.gaze_v_head.r2, regression_stats.head_v_head.r2, regression_stats.pupil_v_head.r2];
tmp_coef_data = [regression_stats.gaze_v_head.coef, regression_stats.head_v_head.coef, regression_stats.pupil_v_head.coef];

figure; set(gcf, 'Position',  [200, 100, 200, 300]); hold on
bar(1:3, mean(tmp_r2_data,1), 'FaceColor', c_head_light)
scatter(1:3, tmp_r2_data, 'k')
plot(1:3, tmp_r2_data, 'k-')
ylim([-0.1 1])
xlim([0.5 3.5])
title('Head regressed on')
ylabel('R2')
xticks([1 2 3])
xticklabels({'Initial gaze', 'Initial head', 'Intial pupil'})
xtickangle(45)

figure; set(gcf, 'Position',  [200, 100, 200, 300]); hold on
bar(1:3, mean(tmp_coef_data,1), 'FaceColor', c_head_light)
scatter(1:3, tmp_coef_data, 'k')
plot(1:3, tmp_coef_data, 'k-')
ylim([-1 1])
xlim([0.5 3.5])
title('Head regressed on')
ylabel('Slope')
xticks([1 2 3])
xticklabels({'Initial gaze', 'Initial head', 'Intial pupil'})
xtickangle(45)


% Pupil amplitude regressed on initial gaze/head/pupil
tmp_r2_data = [regression_stats.gaze_v_pupil.r2, regression_stats.head_v_pupil.r2, regression_stats.pupil_v_pupil.r2];
tmp_coef_data = [regression_stats.gaze_v_pupil.coef, regression_stats.head_v_pupil.coef, regression_stats.pupil_v_pupil.coef];

figure; set(gcf, 'Position',  [200, 100, 200, 300]); hold on
bar(1:3, mean(tmp_r2_data,1), 'FaceColor', c_pupil_light)
scatter(1:3, tmp_r2_data, 'k')
plot(1:3, tmp_r2_data, 'k-')
ylim([-0.1 1])
xlim([0.5 3.5])
title('Pupil regressed on')
ylabel('R2')
xticks([1 2 3])
xticklabels({'Initial gaze', 'Initial head', 'Intial pupil'})
xtickangle(45)

figure; set(gcf, 'Position',  [200, 100, 200, 300]); hold on
bar(1:3, mean(tmp_coef_data,1), 'FaceColor', c_pupil_light)
scatter(1:3, tmp_coef_data, 'k')
plot(1:3, tmp_coef_data, 'k-')
ylim([-1 1])
xlim([0.5 3.5])
title('Pupil regressed on')
ylabel('Slope')
xticks([1 2 3])
xticklabels({'Initial gaze', 'Initial head', 'Intial pupil'})
xtickangle(45)


% Statistics
fprintf('\nR2 values: t test \n')
[~,p] = ttest(regression_stats.gaze_v_gaze.r2, 0);
fprintf('gaze_v_gaze: p: %f\n', p)
[~,p] = ttest(regression_stats.head_v_gaze.r2, 0);
fprintf('head_v_gaze: p: %f\n', p)
[~,p] = ttest(regression_stats.pupil_v_gaze.r2, 0);
fprintf('pupil_v_gaze: p: %f\n', p)

[~,p] = ttest(regression_stats.gaze_v_head.r2, 0);
fprintf('gaze_v_head: p: %f\n', p)
[~,p] = ttest(regression_stats.head_v_head.r2, 0);
fprintf('head_v_head: p: %f\n', p)
[~,p] = ttest(regression_stats.pupil_v_head.r2, 0);
fprintf('pupil_v_head: p: %f\n', p)

[~,p] = ttest(regression_stats.gaze_v_pupil.r2, 0);
fprintf('gaze_v_pupil: p: %f\n', p)
[~,p] = ttest(regression_stats.head_v_pupil.r2, 0);
fprintf('head_v_pupil: p: %f\n', p)
[~,p] = ttest(regression_stats.pupil_v_pupil.r2, 0);
fprintf('pupil_v_pupil: p: %f\n', p)


fprintf('\nSlope values: t test \n')
[~,p] = ttest(regression_stats.gaze_v_gaze.coef, 0);
fprintf('gaze_v_gaze: p: %f\n', p)
[~,p] = ttest(regression_stats.head_v_gaze.coef, 0);
fprintf('head_v_gaze: p: %f\n', p)
[~,p] = ttest(regression_stats.pupil_v_gaze.coef, 0);
fprintf('pupil_v_gaze: p: %f\n', p)

[~,p] = ttest(regression_stats.gaze_v_head.coef, 0);
fprintf('gaze_v_head: p: %f\n', p)
[~,p] = ttest(regression_stats.head_v_head.coef, 0);
fprintf('head_v_head: p: %f\n', p)
[~,p] = ttest(regression_stats.pupil_v_head.coef, 0);
fprintf('pupil_v_head: p: %f\n', p)

[~,p] = ttest(regression_stats.gaze_v_pupil.coef, 0);
fprintf('gaze_v_pupil: p: %f\n', p)
[~,p] = ttest(regression_stats.head_v_pupil.coef, 0);
fprintf('head_v_pupil: p: %f\n', p)
[~,p] = ttest(regression_stats.pupil_v_pupil.coef, 0);
fprintf('pupil_v_pupil: p: %f\n', p)


%% SUPPLEMENT: INTEROCULAR SACCADE CORRELATION

figure; set(gcf, 'Position',  [200, 100, 200, 300]); hold on
bar(1, mean(pupil_saccade_x_correlation), 'FaceColor', c_pupil_light)
scatter(1, pupil_saccade_x_correlation, [], [0, 0, 0])

% Note: Because of our pupil sign convention, a negative correlation in the vertical direction is conjugate movement (e.g., left eye goes down while right eye goes up).
bar(2, mean(pupil_saccade_y_correlation*-1), 'FaceColor', c_pupil_light) % multiplied by -1 so that conjugate movements are positive
scatter(2, pupil_saccade_y_correlation*-1, [], [0, 0, 0]) % multiplied by -1 so that conjugate movements are positive

ylabel('Interocular saccade correlation')
xlim([0, 3])
ylim([-1, 1])
xticks(1:2);
xticklabels({'Horiz.', 'Vert.'});


%% SUPPLEMENT: PLOT MEAN MOVEMENT MAGNITUDES FOR PUPIL HORIZONTAL/VERTICAL

% saccade magnitudes for horizontal and vertical components
figure; set(gcf, 'Position',  [200, 100, 200, 300]); hold on

saccade_horizontal_magnitudes = mean([left_pupil_saccade_amp_mag_x_mean, right_pupil_saccade_amp_mag_x_mean], 2);
saccade_vertical_magnitudes = mean([left_pupil_saccade_amp_mag_y_mean, right_pupil_saccade_amp_mag_y_mean], 2);
bar(1:2, mean([saccade_horizontal_magnitudes, saccade_vertical_magnitudes], 1), 'FaceColor', c_pupil_light)
scatter(1:2, [saccade_horizontal_magnitudes, saccade_vertical_magnitudes], [], [0, 0, 0])

ylabel('Mean saccade magnitude (°)')
xlim([0, 3])
xticks(1:2);
xticklabels({'Horiz.', 'Vert.'});
xtickangle(45)

% statistics
fprintf('\nSaccade magnitude, horizontal: %f +- %f (mean +- SD)\n', mean(saccade_horizontal_magnitudes), std(saccade_horizontal_magnitudes))
fprintf('Saccade magnitude, vertical: %f +- %f (mean +- SD)\n', mean(saccade_vertical_magnitudes), std(saccade_vertical_magnitudes))
[~, p] = ttest2(saccade_horizontal_magnitudes, saccade_vertical_magnitudes);
fprintf('t-test2: p: %f\n', p)


%% SUPPLEMENT: PLOT INITIAL HEAD VS INITIAL PUPIL REGRESSION STATS

figure; set(gcf, 'Position',  [200, 100, 200, 300]); hold on
bar(1, mean(initial_head_vs_initial_pupil.r2, 1), 'FaceColor', [0.5, 0.5, 0.5])
scatter(1, initial_head_vs_initial_pupil.r2, 'k')

ylabel('R2')
xlim([0, 2])

figure; set(gcf, 'Position',  [200, 100, 200, 300]); hold on
bar(1, mean(initial_head_vs_initial_pupil.coef, 1), 'FaceColor', [0.5, 0.5, 0.5])
scatter(1, initial_head_vs_initial_pupil.coef, 'k')

ylabel('Slope')
xlim([0, 2])


%% SUPPLEMENT: PLOT MEAN MOVEMENT MAGNITUDES FOR HEAD YAW/PITCH/ROLL

% head movement magnitudes for yaw, pitch, and roll
figure; set(gcf, 'Position',  [200, 100, 200, 300]); hold on
bar(1:3, [mean(head_yaw_amp_mag_mean), mean(head_roll_amp_mag_mean), mean(head_pitch_amp_mag_mean)], 'FaceColor', c_head_light)
scatter(1:3, [head_yaw_amp_mag_mean, head_roll_amp_mag_mean, head_pitch_amp_mag_mean], [], [0, 0, 0])

ylabel('Mean head mov. magnitude (°)')
xlim([0, 4])
xticks(1:3);
xticklabels({'Yaw', 'Pitch', 'Roll'});
xtickangle(45)

% statistics
fprintf('\nHead yaw magnitude: %f +- %f (mean +- SD)\n', mean(head_yaw_amp_mag_mean), std(head_yaw_amp_mag_mean))
fprintf('Head pitch magnitude: %f +- %f (mean +- SD)\n', mean(head_roll_amp_mag_mean), std(head_roll_amp_mag_mean))
fprintf('Head roll magnitude: %f +- %f (mean +- SD)\n', mean(head_pitch_amp_mag_mean), std(head_pitch_amp_mag_mean))
[~, p] = ttest2(head_yaw_amp_mag_mean, head_roll_amp_mag_mean);
fprintf('t-test2, yaw vs roll: p: %f\n', p)
[~, p] = ttest2(head_yaw_amp_mag_mean, head_pitch_amp_mag_mean);
fprintf('t-test2, yaw vs pitch: p: %f\n', p)
[~, p] = ttest2(head_roll_amp_mag_mean, head_pitch_amp_mag_mean);
fprintf('t-test2, roll vs pitch: p: %f\n', p)


%% SUPPLEMENT: PLOT MAIN SEQUENCE REGRESSION STATS

figure; set(gcf, 'Position',  [200, 100, 200, 300]); hold on
bar(1, mean(saccade_amp_vs_peak_velocity.r2, 1), 'FaceColor', [0.5, 0.5, 0.5])
scatter(1, saccade_amp_vs_peak_velocity.r2, 'k')

ylim([0, 1])
xlim([0, 2])
ylabel('R2')

figure; set(gcf, 'Position',  [200, 100, 200, 300]); hold on
bar(1, mean(saccade_amp_vs_peak_velocity.coef, 1), 'FaceColor', [0.5, 0.5, 0.5])
scatter(1, saccade_amp_vs_peak_velocity.coef, 'k')

xlim([0, 2])
ylabel('Slope')


%% SUPPLEMENT: PLOT AVERAGE INITIAL AND FINAL POSITIONS

tmp_data = [population_averages.initial_head, population_averages.final_head];
figure; set(gcf, 'Position',  [200, 100, 200, 300]); hold on
scatter(1:2, tmp_data, 'k')
plot(1:2, tmp_data', 'k-')
ylim([-40 40])
xlim([0.5 2.5])
ylabel('Position (°)')
xticks([1 2])
xticklabels({'Initial head', 'Final head'})
xtickangle(45)

tmp_data = [population_averages.initial_eye, population_averages.final_eye];
figure; set(gcf, 'Position',  [200, 100, 200, 300]); hold on
scatter(1:2, tmp_data, 'k')
plot(1:2, tmp_data', 'k-')
ylim([-40 40])
xlim([0.5 2.5])
ylabel('Position (°)')
xticks([1 2])
xticklabels({'Initial pupil', 'Final pupil'})
xtickangle(45)


%% SUPPLEMENT: PLOT AVERAGE AMPLITUDES

tmp_data = population_averages.amplitude_head;
figure; set(gcf, 'Position',  [200, 100, 100, 300]); hold on
scatter(1, tmp_data, 'k')
ylim([-40 40])
xlim([0.5 1.5])
ylabel('Mean head displacement (°)')
xticks([1])

tmp_data = population_averages.amplitude_eye;
figure; set(gcf, 'Position',  [200, 100, 100, 300]); hold on
scatter(1, tmp_data, 'k')
ylim([-40 40])
xlim([0.5 1.5])
ylabel('Mean saccade displacement (°)')
xticks([1])

