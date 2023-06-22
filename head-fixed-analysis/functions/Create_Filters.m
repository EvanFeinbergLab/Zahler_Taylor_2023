function [Filters,hankelMtrx] = Create_Filters(mtrx,trial_param,filter_param)
%Create_Filters takes a matrix of time series data and separates it into
%trials and outputs a hankelMtrx as well as logical descriptors of each
%trial in Filters

if isfield(filter_param, 'pupil_side')
    pupil_side = filter_param.pupil_side;
else
    pupil_side = 2;
end

if isfield(filter_param, 'strain_gauge_movement_threshold')
    strain_gauge_movement_threshold = filter_param.strain_gauge_movement_threshold;
else
    strain_gauge_movement_threshold = 0.2;
end

if isfield(filter_param,'smooth_strain_gauge_window')
    smooth_strain_gauge_window = filter_param.smooth_strain_gauge_window;
else
    smooth_strain_gauge_window = 1:trial_param.align-1; 
end

if isfield(filter_param,'strain_gauge_after_align_window')
    strain_gauge_after_align_window = filter_param.strain_gauge_after_align_window;
else
    strain_gauge_after_align_window = trial_param.align:trial_param.align+8; 
end

if isfield(filter_param,'small_pupil_window')
    small_pupil_window = filter_param.small_pupil_window;
else
    small_pupil_window = 1:trial_param.align-1;
end

if isfield(filter_param,'small_pupil_threshold')
    small_pupil_threshold = filter_param.small_pupil_threshold;
else
    small_pupil_threshold = 0;
end

if isfield(filter_param,'no_saccade_before_align_period')
    no_saccade_before_align_period = filter_param.no_saccade_before_align_period;
else
    no_saccade_before_align_period = trial_param.preAlignPeriod;
end

if isfield(filter_param, 'pupil_movement_threshold')
    pupil_movement_threshold = filter_param.pupil_movement_threshold;
else
    pupil_movement_threshold = 5;
end

if isfield(filter_param,'smooth_pupil_window')
    smooth_pupil_window = filter_param.smooth_pupil_window;
else
    smooth_pupil_window = 1:trial_param.align-1; 
end

%Identify and align trial starts
starts_shift = trial_param.starts-trial_param.align+1;

% mouse id for each trial
tmp_mouseID = align2starts(mtrx.mouseID, starts_shift,trial_param.window_size);
mouseID = tmp_mouseID(:,trial_param.align);

%Matrix of aligned pupil traces. Use L/R/Mean pupil depending on pupil_side
if pupil_side == 0
    [unalignedPupil, starts_shift_trim] = align2starts(mtrx.pupilX_L, starts_shift,trial_param.window_size);
elseif pupil_side == 1
    [unalignedPupil, starts_shift_trim] = align2starts(mtrx.pupilX_R, starts_shift,trial_param.window_size);
elseif pupil_side == 2
    [unalignedPupil, starts_shift_trim] = align2starts(mtrx.pupilX_Mean, starts_shift,trial_param.window_size);
end

% pupil elevation
if pupil_side == 0
    [unalignedPupilY, starts_shift_trim] = align2starts(mtrx.pupilY_L, starts_shift,trial_param.window_size);
elseif pupil_side == 1
    [unalignedPupilY, starts_shift_trim] = align2starts(mtrx.pupilY_R, starts_shift,trial_param.window_size);
elseif pupil_side == 2
    [unalignedPupilY, starts_shift_trim] = align2starts(mtrx.pupilY_Mean, starts_shift,trial_param.window_size);
end

%pupil position at align
initPupilPos = unalignedPupil(:,trial_param.align);
%pupil position at end of response window
finalPupilPos = unalignedPupil(:,trial_param.response_window(length(trial_param.response_window)));
alignedPupil = unalignedPupil-initPupilPos;

%pupil position at align (elevation)
initPupilPosY = unalignedPupilY(:,trial_param.align);
%pupil position at end of response window (elevation)
finalPupilPosY = unalignedPupilY(:,trial_param.response_window(length(trial_param.response_window)));
alignedPupilY = unalignedPupilY-initPupilPosY;

% matrix of saccade amplitudes at saccade onset
saccade = align2starts(mtrx.saccade,starts_shift_trim,trial_param.window_size);

% matrix of saccade amplitudes at saccade onset
saccade_length = align2starts(mtrx.saccade_length,starts_shift_trim,trial_param.window_size);

% find_saccades_2d
saccade_2d_angle = align2starts(mtrx.saccade_2d_angle,starts_shift_trim,trial_param.window_size);
saccade_2d_amp = align2starts(mtrx.saccade_2d_amp,starts_shift_trim,trial_param.window_size);
saccade_2d_x = align2starts(mtrx.saccade_2d_x,starts_shift_trim,trial_param.window_size);
saccade_2d_y = align2starts(mtrx.saccade_2d_y,starts_shift_trim,trial_param.window_size);

% matrix of strain gauge traces 
unalignedStrainGauge = align2starts(mtrx.strain_gauge,starts_shift_trim,trial_param.window_size);

if isfield(trial_param,'alignStrainGauge')
    alignStrainGauge = trial_param.alignStrainGauge;
else
    alignStrainGauge = 1; 
end

initStrainGaugePos = unalignedStrainGauge(:,alignStrainGauge);
alignedStrainGauge = unalignedStrainGauge-initStrainGaugePos;

% NULL FILTER
numStarts = size(unalignedPupil,1);
Filters.null = ones(numStarts,1,'logical');

% LEFT/RIGHT PUFF FILTERS
leftPuff = align2starts(mtrx.leftPuff,starts_shift_trim,trial_param.window_size);
Filters.leftPuffFilter = sum(leftPuff,2);
rightPuff = align2starts(mtrx.rightPuff,starts_shift_trim,trial_param.window_size);
Filters.rightPuffFilter = sum(rightPuff,2);

% LEFT2/RIGHT2 PUFF FILTERS
leftPuff2 = align2starts(mtrx.leftPuff2,starts_shift_trim,trial_param.window_size);
Filters.leftPuff2Filter = sum(leftPuff2,2);
rightPuff2 = align2starts(mtrx.rightPuff2,starts_shift_trim,trial_param.window_size);
Filters.rightPuff2Filter = sum(rightPuff2,2);

% ALL PUFF FILTERS
Filters.allPuffFilter = Filters.leftPuffFilter | Filters.rightPuffFilter | Filters.leftPuff2Filter | Filters.rightPuff2Filter;


% SMOOTH STRAIN GAUGE FILTER
% Strain gauge activity must be constant during a specified window of time.
numStarts = size(unalignedPupil,1);
Filters.smoothStrainGaugeFilter = zeros(numStarts,1,'logical');
if isfield(mtrx,'strain_gauge')
    for i = 1:size(Filters.smoothStrainGaugeFilter,1)
        if isempty(find(range(unalignedStrainGauge(i,smooth_strain_gauge_window))>strain_gauge_movement_threshold, 1))
            Filters.smoothStrainGaugeFilter(i) = 1;
        end   
    end
end

% STRAIN GAUGE AFTER ALIGN FILTER
numStarts = size(unalignedPupil,1);
Filters.strainGaugeAfterAlignFilter = zeros(numStarts,1,'logical');
if isfield(mtrx,'strain_gauge')
    for i = 1:size(Filters.strainGaugeAfterAlignFilter,1)
        if ~isempty(find(range(unalignedStrainGauge(i,strain_gauge_after_align_window))>strain_gauge_movement_threshold, 1))
            Filters.strainGaugeAfterAlignFilter(i) = 1;
        end   
    end
end

%SMALL PUPIL FILTER
%Pupil diameter must be smaller than threshold value
Filters.smallPupilFilter = zeros(numStarts,1,'logical');
pupilArea = align2starts(mtrx.pupilArea,starts_shift,trial_param.window_size);
for i = 1:size(Filters.smallPupilFilter,1)
    if isempty(find(pupilArea(i,small_pupil_window)>small_pupil_threshold, 1))
        Filters.smallPupilFilter(i) = 1;
    end
end

% SMOOTH PUPIL FILTER
Filters.smoothPupilFilter = ~any(range(unalignedPupil(:,smooth_pupil_window), 2)>pupil_movement_threshold, 2);

% SACCADE BEFORE ALIGN
Filters.saccadeBeforeAlignFilter = zeros(numStarts,1,'logical');
for i = 1:length(Filters.saccadeBeforeAlignFilter)
    if ~isempty(find(saccade(i,no_saccade_before_align_period)~=0, 1))
        Filters.saccadeBeforeAlignFilter(i) = 1;
    end
end

%SACCADE DURING RESPONSE WINDOW FILTER
Filters.saccadeAfterAlignFilter = zeros(numStarts,1,'logical');
Filters.leftSaccadeAfterAlignFilter = zeros(numStarts,1,'logical');
Filters.rightSaccadeAfterAlignFilter = zeros(numStarts,1,'logical');
for i = 1:numStarts
    if ~isempty(find(saccade(i,trial_param.response_window)~=0, 1))
        Filters.saccadeAfterAlignFilter(i) = 1; 
    end 
    if ~isempty(find(saccade(i,trial_param.response_window)<0, 1))
        Filters.leftSaccadeAfterAlignFilter(i) = 1;
    end 
    if ~isempty(find(saccade(i,trial_param.response_window)>0, 1))
        Filters.rightSaccadeAfterAlignFilter(i) = 1;
    end
end

% SACCADE AFTER RESPONSE WINDOW FILTER
Filters.lateSaccadeAfterAlignFilter = zeros(numStarts,1,'logical');
for i = 1:numStarts
    if ~isempty(find(saccade(i,trial_param.response_window(end)+1:end)~=0, 1))
        Filters.lateSaccadeAfterAlignFilter(i) = 1;
    end 
end

% NaN DURING TRIAL FILTER
% NaNs occur when dlc tracking likelihood is below 0.99 for any point
Filters.noNanFilter = zeros(numStarts,1,'logical');
for i = 1:numStarts
    if ~ismember(1,isnan(unalignedPupil(i,:))) && ~ismember(1,isnan(unalignedStrainGauge(i,:)))
        Filters.noNanFilter(i) = 1;
    end
end

Filters.includeAllFilter = ones(numStarts,1,'logical');

%OPTO STIM DURING TRIAL FILTER
Filters.optoTrialFilter = zeros(numStarts,1,'logical');
tmpOpto = align2starts(mtrx.optoTrial,starts_shift_trim,trial_param.window_size);
for i = 1:numStarts
    if ~isempty(find(tmpOpto(i,:) == 1, 1))
        Filters.optoTrialFilter(i) = 1;
    end
end

if isfield(filter_param, 'stable_recording_period')
    Filters.stable_recording_period = starts_shift_trim < (filter_param.stable_recording_period*trial_param.frame_rate);
end

%OPTO INTENSITY (OG BOX UNITS) DURING TRIAL FILTER
tmpOptoIntensity = align2starts(mtrx.opto_intensity,starts_shift_trim,trial_param.window_size);
Filters.optoIntensityFilter = sum(tmpOptoIntensity,2);

%GENERATE OUTPUT
hankelMtrx.mouseID = mouseID; % encoded mouse id for each trial
hankelMtrx.saccade = saccade;
hankelMtrx.saccade_length = saccade_length;
hankelMtrx.unalignedPupil = unalignedPupil; 
hankelMtrx.alignedPupil = alignedPupil; 
hankelMtrx.unalignedStrainGauge = unalignedStrainGauge;
hankelMtrx.alignedStrainGauge = alignedStrainGauge;
hankelMtrx.initPupilPos = initPupilPos;
hankelMtrx.finalPupilPos = finalPupilPos;
hankelMtrx.pupilArea = pupilArea;

% find_saccades_2d
hankelMtrx.saccade_2d_angle = saccade_2d_angle;
hankelMtrx.saccade_2d_amp = saccade_2d_amp;
hankelMtrx.saccade_2d_x = saccade_2d_x;
hankelMtrx.saccade_2d_y = saccade_2d_y;

% pupilY
hankelMtrx.unalignedPupilY = unalignedPupilY; 
hankelMtrx.alignedPupilY = alignedPupilY; 
hankelMtrx.initPupilPosY = initPupilPosY;
hankelMtrx.finalPupilPosY = finalPupilPosY;

% puff 
hankelMtrx.leftPuff = Filters.leftPuffFilter;
hankelMtrx.rightPuff = Filters.rightPuffFilter;


% opto
hankelMtrx.optoOnset = align2starts(mtrx.optoTrial, starts_shift,trial_param.window_size);

end

