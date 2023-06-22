function [Filters,hankelMtrx] = Create_Filters(mtrx,trial_param,filter_param)
%Create_Filters takes a matrix of time series data and separates it into
%trials and outputs a hankelMtrx as well as logical descriptors of each
%trial in Filters

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

if isfield(filter_param, 'head_movement_threshold')
    head_movement_threshold = filter_param.head_movement_threshold;
else
    head_movement_threshold = 0.2;
end

if isfield(filter_param,'smooth_head_window')
    smooth_head_window = filter_param.smooth_head_window;
else
    smooth_head_window = 1:trial_param.align-1; 
end

%Identify and align trial starts
starts_shift = trial_param.starts-trial_param.align+1;

% CREATE HANKEL MATRICES

% mouse id for each trial
tmp_mouseID = align2starts(mtrx.mouseID, starts_shift,trial_param.window_size);
mouseID = tmp_mouseID(:,trial_param.align);

hankelMtrx.mouseID = mouseID; % encoded mouse id for each trial

% Opto onset hankel matrix
hankelMtrx.optoOnset = align2starts(mtrx.optoTrial, starts_shift,trial_param.window_size);

% Airpuff onset hankel matrix
hankelMtrx.leftAirpuffOnset = align2starts(mtrx.leftPuff, starts_shift,trial_param.window_size);
hankelMtrx.rightAirpuffOnset = align2starts(mtrx.rightPuff, starts_shift,trial_param.window_size);

% PUPIL HANKEL MATRICES AND FILTERS
if ~isempty(mtrx.pupilX_L)
    hankelMtrx.leftPupilX = align2starts(mtrx.pupilX_L, starts_shift,trial_param.window_size);
    hankelMtrx.leftPupilY = align2starts(mtrx.pupilY_L, starts_shift,trial_param.window_size);
    
    hankelMtrx.rightPupilX = align2starts(mtrx.pupilX_R, starts_shift,trial_param.window_size);
    hankelMtrx.rightPupilY = align2starts(mtrx.pupilY_R, starts_shift,trial_param.window_size);
    
    hankelMtrx.meanPupilX = align2starts(mtrx.pupilX_Mean, starts_shift,trial_param.window_size);
    hankelMtrx.meanPupilY = align2starts(mtrx.pupilY_Mean, starts_shift,trial_param.window_size);



    %SMALL PUPIL FILTER
    %Pupil diameter must be smaller than threshold value at all timepoints
    %within small_pupil_window
    hankelMtrx.pupilArea = align2starts(mtrx.pupilArea,starts_shift,trial_param.window_size);
    Filters.smallPupilFilter = ~any(hankelMtrx.pupilArea(:,small_pupil_window)>small_pupil_threshold, 2);
    
    % NaN DURING TRIAL FILTER
    % NaNs occur when dlc tracking likelihood falls below threshold
    % (threshold set in combine_sessions_freely_moving function)
    Filters.noLeftPupilNaNFilter = ~(any(isnan(hankelMtrx.leftPupilX), 2));
    Filters.noRightPupilNaNFilter = ~(any(isnan(hankelMtrx.rightPupilX), 2));
    
    % SACCADES HANKEL
    hankelMtrx.leftPupilXSaccades = align2starts(mtrx.pupilX_L_Saccade, starts_shift,trial_param.window_size);
    hankelMtrx.leftPupilYSaccades = align2starts(mtrx.pupilY_L_Saccade, starts_shift,trial_param.window_size);
    
    hankelMtrx.rightPupilXSaccades = align2starts(mtrx.pupilX_R_Saccade, starts_shift,trial_param.window_size);
    hankelMtrx.rightPupilYSaccades = align2starts(mtrx.pupilY_R_Saccade, starts_shift,trial_param.window_size);
    
    hankelMtrx.leftPupilSaccadesLength = align2starts(mtrx.pupil_L_Saccade_Length, starts_shift,trial_param.window_size);
    hankelMtrx.rightPupilSaccadesLength = align2starts(mtrx.pupil_R_Saccade_Length, starts_shift,trial_param.window_size);

    hankelMtrx.leftPupilSaccadesPeakVelocityOffset = align2starts(mtrx.pupilX_L_Saccade_PeakVelocityOffset, starts_shift,trial_param.window_size);
    hankelMtrx.rightPupilSaccadesPeakVelocityOffset = align2starts(mtrx.pupilX_R_Saccade_PeakVelocityOffset, starts_shift,trial_param.window_size);

    hankelMtrx.leftPupilSaccadesPeakVelocity = align2starts(mtrx.pupilX_L_Saccade_PeakVelocity, starts_shift,trial_param.window_size);
    hankelMtrx.rightPupilSaccadesPeakVelocity = align2starts(mtrx.pupilX_R_Saccade_PeakVelocity, starts_shift,trial_param.window_size);
    
    %SACCADE DURING RESPONSE WINDOW FILTER
    Filters.leftPupilSaccadeAfterAlignFilter = any(hankelMtrx.leftPupilXSaccades(:,trial_param.response_window),2);
    Filters.leftPupilLeftSaccadeAfterAlignFilter = Filters.leftPupilSaccadeAfterAlignFilter & row_indexing(hankelMtrx.leftPupilXSaccades, first_event_after(hankelMtrx.leftPupilXSaccades, trial_param.align))<0;
    Filters.leftPupilRightSaccadeAfterAlignFilter = Filters.leftPupilSaccadeAfterAlignFilter & row_indexing(hankelMtrx.leftPupilXSaccades, first_event_after(hankelMtrx.leftPupilXSaccades, trial_param.align))>0;
    
    Filters.rightPupilSaccadeAfterAlignFilter = any(hankelMtrx.rightPupilXSaccades(:,trial_param.response_window),2);
    Filters.rightPupilLeftSaccadeAfterAlignFilter = Filters.rightPupilSaccadeAfterAlignFilter & row_indexing(hankelMtrx.rightPupilXSaccades, first_event_after(hankelMtrx.rightPupilXSaccades, trial_param.align))<0;
    Filters.rightPupilRightSaccadeAfterAlignFilter = Filters.rightPupilSaccadeAfterAlignFilter & row_indexing(hankelMtrx.rightPupilXSaccades, first_event_after(hankelMtrx.rightPupilXSaccades, trial_param.align))>0;
    
    Filters.bothPupilSaccadeAfterAlignFilter = Filters.leftPupilSaccadeAfterAlignFilter & Filters.rightPupilSaccadeAfterAlignFilter;
    Filters.bothPupilLeftSaccadeAfterAlignFilter = Filters.leftPupilLeftSaccadeAfterAlignFilter & Filters.rightPupilLeftSaccadeAfterAlignFilter;
    Filters.bothPupilRightSaccadeAfterAlignFilter = Filters.leftPupilRightSaccadeAfterAlignFilter & Filters.rightPupilRightSaccadeAfterAlignFilter;
    
    % simultaneousSaccadeAfterAlignFilter
    Filters.simultaneousSaccadeAfterAlignFilter = logical(zeros(size(hankelMtrx.leftPupilXSaccades, 1), 1));
    for i = 1:size(hankelMtrx.leftPupilXSaccades, 1)
        tmp_left_pupil_first_saccade_after_idx = first_event_after(hankelMtrx.leftPupilXSaccades(i,:), trial_param.align);
        tmp_right_pupil_first_saccade_after_idx = first_event_after(hankelMtrx.rightPupilXSaccades(i,:), trial_param.align);
        
        if any(trial_param.response_window == tmp_left_pupil_first_saccade_after_idx) && any(trial_param.response_window == tmp_right_pupil_first_saccade_after_idx) % must have saccades in response window
            if abs(tmp_left_pupil_first_saccade_after_idx - tmp_right_pupil_first_saccade_after_idx) <= 3  % saccades must be within 3 frames of each other
                if sign(hankelMtrx.leftPupilXSaccades(i,tmp_left_pupil_first_saccade_after_idx)) == sign(hankelMtrx.rightPupilXSaccades(i,tmp_right_pupil_first_saccade_after_idx)) % saccades must be in same direction
                    Filters.simultaneousSaccadeAfterAlignFilter(i) = 1;
                end
            end
        end
    end


    % SMOOTH PUPIL FILTERS
    % These filters identify trials where the pupil trace stays relatively constant before saccade onset. 
    % The range of pupil trace values within a window preceding trial onset (i.e., smooth_pupil_window)
    % cannot exceed pupil_movement_threshold.
    Filters.smoothLeftPupilXFilter = ~any(range(hankelMtrx.leftPupilX(:,smooth_pupil_window), 2)>pupil_movement_threshold, 2);
    Filters.smoothLeftPupilYFilter = ~any(range(hankelMtrx.leftPupilY(:,smooth_pupil_window), 2)>pupil_movement_threshold, 2);
    Filters.smoothRightPupilXFilter = ~any(range(hankelMtrx.rightPupilX(:,smooth_pupil_window), 2)>pupil_movement_threshold, 2);
    Filters.smoothRightPupilYFilter = ~any(range(hankelMtrx.rightPupilY(:,smooth_pupil_window), 2)>pupil_movement_threshold, 2);
    Filters.smoothLeftPupilFilter = Filters.smoothLeftPupilXFilter & Filters.smoothLeftPupilYFilter;
    Filters.smoothRightPupilFilter = Filters.smoothRightPupilXFilter & Filters.smoothRightPupilYFilter;
    Filters.smoothPupilFilter = Filters.smoothLeftPupilFilter & Filters.smoothRightPupilFilter; % this filter rejects pre-trial movements of left and right pupils in both vertical and horizontal axes


    % NO MOVEMENT BEFORE FIRST SACCADE FILTER
    % This filter identifies trials where the difference between the pupil value at trial start (align) and saccade onset is less than
    % pupil_movement_threshold. This is intended to eliminate trials where there is a VOR movement before a saccade
    num_trials = size(hankelMtrx.leftPupilXSaccades, 1);
    
    % initialize arrays to store size of pre-saccade movement
    tmp_left_pupil_x_presaccade_movement = NaN(num_trials, 1);
    tmp_right_pupil_x_presaccade_movement = NaN(num_trials, 1);
    
    % Loop through each trial. Record pre-saccade movement for trials that have saccades
    for TRIAL = 1:num_trials
        tmp_first_saccade_after_idx = first_event_after(hankelMtrx.leftPupilXSaccades(TRIAL,:), trial_param.align);
        if ~isnan(tmp_first_saccade_after_idx)
            tmp_left_pupil_x_presaccade_movement(TRIAL,1) = hankelMtrx.leftPupilX(TRIAL,trial_param.align) - hankelMtrx.leftPupilX(TRIAL,tmp_first_saccade_after_idx);
        end
        
        tmp_first_saccade_after_idx = first_event_after(hankelMtrx.rightPupilXSaccades(TRIAL,:), trial_param.align);
        if ~isnan(tmp_first_saccade_after_idx)
            tmp_right_pupil_x_presaccade_movement(TRIAL,1) = hankelMtrx.rightPupilX(TRIAL,trial_param.align) - hankelMtrx.rightPupilX(TRIAL,tmp_first_saccade_after_idx);
        end
    end
    Filters.noMovementBeforeFirstSaccadeFilter = all(abs([tmp_left_pupil_x_presaccade_movement, tmp_right_pupil_x_presaccade_movement]) < pupil_movement_threshold, 2); % both pupils must fall below the threshold

end


% IMU HANKEL MATRICES AND FILTERS
if isfield(mtrx, 'rawImuData')
    

    % Add raw IMU data to hankelMtrx
    for i = 1:numel(fieldnames(mtrx.rawImuData))
        mtrx_fieldnames = fieldnames(mtrx.rawImuData);
        rawImuData.(mtrx_fieldnames{i}) = align2starts(mtrx.rawImuData.(mtrx_fieldnames{i}),starts_shift,trial_param.window_size);
    end    
    hankelMtrx.rawImuData = rawImuData;
    
    % Assign processed IMU data to hankelMtrx
    hankelMtrx.headWorldYawVelIMU = align2starts(mtrx.headWorldYawVelIMU,starts_shift,trial_param.window_size);
    hankelMtrx.headWorldYawAngleIMU = align2starts(mtrx.headWorldYawAngleIMU,starts_shift,trial_param.window_size); % use session velocity cumsum
%     hankelMtrx.headWorldYawAngleIMU = cumsum(hankelMtrx.headWorldYawVelIMU, 2); % alternative: cumsum velocity each trial
    hankelMtrx.headPitchAngle = align2starts(mtrx.headPitchAngle,starts_shift,trial_param.window_size);
    hankelMtrx.headRollAngle = align2starts(mtrx.headRollAngle,starts_shift,trial_param.window_size);

    hankelMtrx.headYawEventsAmp = align2starts(mtrx.headYawEvents.amp,starts_shift,trial_param.window_size);
    hankelMtrx.headYawEventsDuration = align2starts(mtrx.headYawEvents.duration,starts_shift,trial_param.window_size);

    hankelMtrx.headWorldPitchVelGyro = align2starts(mtrx.headWorldPitchVelGyro,starts_shift,trial_param.window_size);
    hankelMtrx.headWorldPitchAngleGyro = align2starts(mtrx.headWorldPitchAngleGyro,starts_shift,trial_param.window_size); % use session velocity cumsum
    
    hankelMtrx.headWorldRollVelGyro = align2starts(mtrx.headWorldRollVelGyro,starts_shift,trial_param.window_size);
    hankelMtrx.headWorldRollAngleGyro = align2starts(mtrx.headWorldRollAngleGyro,starts_shift,trial_param.window_size); % use session velocity cumsum

    % IMU FILTERS

    % SMOOTH HEAD FILTER (YAW)
    % The range of head trace values within a window preceding trial onset (i.e., smooth_head_window) cannot exceed head_movement_threshold.
    Filters.smoothHeadFilter = ~any(range(hankelMtrx.headWorldYawVelIMU(:,smooth_head_window), 2)>head_movement_threshold, 2);
    
    % HEAD MOVEMENT EVENT FILTERS
    Filters.headYawEventAfterAlignFilter = any(hankelMtrx.headYawEventsAmp(:,trial_param.response_window),2); % true if any yaw head movement event found in response window
    Filters.leftHeadYawEventAfterAlignFilter = Filters.headYawEventAfterAlignFilter & row_indexing(hankelMtrx.headYawEventsAmp, first_event_after(hankelMtrx.headYawEventsAmp, trial_param.align))<0; % true if any left yaw head movement event found in response window
    Filters.rightHeadYawEventAfterAlignFilter = Filters.headYawEventAfterAlignFilter & row_indexing(hankelMtrx.headYawEventsAmp, first_event_after(hankelMtrx.headYawEventsAmp, trial_param.align))>0; % true if any right yaw head movement event found in response window
    
end

% BODY CAMERA HANKEL MATRICES AND FILTERS
if isfield(mtrx, 'rawBodyCamData')
    for i = 1:numel(fieldnames(mtrx.rawBodyCamData))
        mtrx_fieldnames = fieldnames(mtrx.rawBodyCamData);
        rawBodyCamData.(mtrx_fieldnames{i}) = align2starts(mtrx.rawBodyCamData.(mtrx_fieldnames{i}),starts_shift,trial_param.window_size);
    end    

    hankelMtrx.rawBodyCamData = rawBodyCamData;
    
    hankelMtrx.headBodyYawAngleVid = align2starts(mtrx.headBodyYawAngleVid,starts_shift,trial_param.window_size);
    hankelMtrx.bodyWorldYawAngleVid = align2starts(mtrx.bodyWorldYawAngleVid,starts_shift,trial_param.window_size);
    hankelMtrx.platformAngle = align2starts(mtrx.platformAngle,starts_shift,trial_param.window_size);
    
    % BODYCAM NaN DURING TRIAL FILTER
    Filters.noBodycamNaNFilter = ~(any(isnan(hankelMtrx.headBodyYawAngleVid), 2));
    
    % NO PLATFORM MOVEMENT FILTER
    Filters.noPlatformMovementFilter = all(abs(hankelMtrx.platformAngle) > 165, 2) | all(abs(hankelMtrx.platformAngle) < 15, 2);
end


% Set headBodyYawAngleIMU starting positions using headBodyYawAngleVid
if isfield(mtrx, 'rawBodyCamData') && isfield(mtrx, 'rawImuData')
    hankelMtrx.headBodyYawAngleIMU = hankelMtrx.headWorldYawAngleIMU - hankelMtrx.headWorldYawAngleIMU(:,trial_param.align) + hankelMtrx.headBodyYawAngleVid(:,trial_param.align);

    hankelMtrx.headBodyPitchAngleGyro = hankelMtrx.headWorldPitchAngleGyro - hankelMtrx.headWorldPitchAngleGyro(:,trial_param.align) + hankelMtrx.headPitchAngle(:,trial_param.align);
    hankelMtrx.headBodyRollAngleGyro = hankelMtrx.headWorldRollAngleGyro - hankelMtrx.headWorldRollAngleGyro(:,trial_param.align) + hankelMtrx.headRollAngle(:,trial_param.align);
end


% CREATE ADDITIONAL FILTERS

% NULL FILTER
numStarts = numel(starts_shift);
Filters.includeAllFilter = ones(numStarts,1,'logical');

% LEFT/RIGHT PUFF FILTERS
leftPuff = align2starts(mtrx.leftPuff,starts_shift,trial_param.window_size);
Filters.leftPuffFilter = sum(leftPuff,2);
rightPuff = align2starts(mtrx.rightPuff,starts_shift,trial_param.window_size);
Filters.rightPuffFilter = sum(rightPuff,2);

% ALL PUFF FILTERS
Filters.allPuffFilter = Filters.leftPuffFilter | Filters.rightPuffFilter;

% NaN DURING TRIAL FILTER
% NaNs occur when dlc tracking likelihood falls below threshold
% (threshold set in combine_sessions_freely_moving function)

if isfield(Filters, 'noBodycamNaNFilter') && isfield(Filters, 'noLeftPupilNaNFilter') && isfield(Filters, 'noRightPupilNaNFilter')
    Filters.noNanFilter = Filters.noBodycamNaNFilter & Filters.noLeftPupilNaNFilter & Filters.noRightPupilNaNFilter;
%     Filters.noNanFilter = Filters.noLeftPupilNaNFilter & Filters.noRightPupilNaNFilter;
elseif isfield(Filters, 'noBodycamNaNFilter')
    Filters.noNanFilter = Filters.noBodycamNaNFilter;
elseif isfield(Filters, 'noLeftPupilNaNFilter') && isfield(Filters, 'noRightPupilNaNFilter')
    Filters.noNanFilter = Filters.noLeftPupilNaNFilter & Filters.noRightPupilNaNFilter;
else
    Filters.noNanFilter = Filters.includeAllFilter;
end

%OPTO STIM DURING TRIAL FILTER
tmpOpto = align2starts(mtrx.optoTrial,starts_shift,trial_param.window_size);
Filters.optoTrialFilter = any(tmpOpto, 2);

%OPTO INTENSITY (OG BOX UNITS) DURING TRIAL FILTER
tmpOptoIntensity = align2starts(mtrx.opto_intensity,starts_shift,trial_param.window_size);
Filters.optoIntensityFilter = sum(tmpOptoIntensity,2);



end

