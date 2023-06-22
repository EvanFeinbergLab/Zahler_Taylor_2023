function [ saccade, saccade_length, velocity, saccade_seed ] = find_saccades_freely_moving(pupil_position, saccade_amplitude_threshold, saccade_velocity_threshold, saccade_onset_velocity_threshold, saccade_offset_velocity_threshold, frame_rate)
% find_saccades_freely_moving identifies the start, end, and amplitude of saccades, optomized for freely moving data
% Adapted from Meyer et al., 2020

% Overview:
% Seed positions for tentative saccades are identified as timepoints where pupil velocity exceeded saccade_velocity_threshold (°/s) in the horizontal axis. 
% To ensure that saccades are not counted more than once, a seed timepoint is only considered if it is the local maximum velocity within a ±50 ms time window. 
% To identify the saccade onset and offset, a ±50 ms window around each seed is oversampled by a factor of two (5 ms resolution). 
% Saccade onset and offset are defined as the first timepoint where one of the following criteria are met: velocity falls below 
% saccade_onset_velocity_threshold/saccade_offset_velocity_threshold (°/s), velocity changes direction, or pupil movement changes direction.
% If none of those criteria are met, onset or offset is set at ±35 ms.

% INPUTS
% 
% pupil_position: Mx1 array of pupil position (degrees)
% saccade_amplitude_threshold: NOT IN USE
% saccade_velocity_threshold: (deg/sec)
% saccade_onset_velocity_threshold:  (deg/sec)
% saccade_offset_velocity_threshold:  (deg/sec)
% frame_rate: frames per second
%
% OUTPUTS
%
% saccade: Mx1 sparse array of saccade amplitudes. Values located at saccade onset frames.
% saccade_length: Mx1 sparse array of saccade lengths. Values located at saccade onset frames.
% velocity: Mx1 array of pupil velocity (deg/sec)
% saccade_seed: Mx1 sparse array of seed location (i.g., peak velocity) relative to saccade onset. Values located at saccade onset frames.


% invert pupil_position if M<N
if size(pupil_position,1) < size(pupil_position,2)
    pupil_position = pupil_position';
end

%calculate time between frames in seconds
dt = 1/frame_rate;

% threshold velocity to find tentative saccades
velocity = [0; diff(pupil_position)]/dt; 
velocityBool = abs(velocity)>saccade_velocity_threshold; %identify frames in which velocity exceeds velocity threshold

%initialize saccade arrays
saccade = zeros(length(velocity),1);
saccade_length = zeros(length(velocity),1);
saccade_seed = zeros(length(velocity),1);

% define window to look for local maxima
local_extrema_window_sec = 0.05; % Seconds. Default = 0.05. For each tentative saccade found by thresholding, look for local peak within this window. 
local_extrema_window_frames = local_extrema_window_sec/dt; % convert to frames

for CURRENT_FRAME = local_extrema_window_frames+1:length(velocityBool)-local_extrema_window_frames-1
    
    % Check if velocity at CURRENT_FRAME exceeds velocity threshold
    if velocityBool(CURRENT_FRAME) == 1 
        
        % Check if velocity at CURRENT_FRAME is local maximum (using local_extrema_window_size)
        tmp_local_vel_abs = abs(velocity(CURRENT_FRAME-local_extrema_window_frames:CURRENT_FRAME+local_extrema_window_frames));
        [~, local_max_idx] = max(tmp_local_vel_abs);
        if local_max_idx == local_extrema_window_frames + 1
            
            % oversample velocity trace through interpolation
            oversampling_factor = 2; % larger values increase time resolution (e.g., if normal time resolution is 10ms, oversampling with a factor of 2 increases resolution to 5ms)
            dt_oversampled = dt/oversampling_factor;
            
            v = pupil_position(CURRENT_FRAME-local_extrema_window_frames:CURRENT_FRAME+local_extrema_window_frames);
            x = 1:numel(v);
            xq = linspace(1,1+2*local_extrema_window_frames,local_extrema_window_frames*2*oversampling_factor+1);
            vq = interp1(x,v,xq);
            vq_vel = diff_pad(vq,2)/dt_oversampled;
            middle = round(numel(vq_vel)/2);
            
            onset_boundary = 0.035/dt_oversampled; % if none of the onset/offset criteria are met, the algorithm stops at +- 35 ms

            % find onset
            for j = -1:-1:-onset_boundary
                velocity_below_threshold = abs(vq_vel(middle+j)) < saccade_onset_velocity_threshold; % criteria 1: velocity falls below threshold
                velocity_sign_change = sign(vq_vel(middle+j-1)) ~= sign(vq_vel(middle+j)); % criteria 2: velocity changes direction
                direction_change = sign(diff(vq(middle+j-1:middle+j)-vq(middle+j))) ~= sign(diff(vq(middle+j:middle+j+1)-vq(middle+j))); % criteria 2: pupil movement changes direction
                if velocity_below_threshold || velocity_sign_change || direction_change
                    onset_idx = CURRENT_FRAME+round(j/oversampling_factor);
                    break;
                end
                onset_idx = CURRENT_FRAME+round(j/oversampling_factor);
            end
            
            % find offset
            for j = 1:1:onset_boundary
                velocity_below_threshold = abs(vq_vel(middle+j)) < saccade_offset_velocity_threshold; % criteria 1: velocity falls below threshold
                velocity_sign_change = sign(vq_vel(middle+j-1)) ~= sign(vq_vel(middle+j)); % criteria 2: velocity changes direction
                direction_change = sign(diff(vq(middle+j-1:middle+j)-vq(middle+j))) ~= sign(diff(vq(middle+j:middle+j+1)-vq(middle+j))); % criteria 2: pupil movement changes direction
                if velocity_below_threshold || velocity_sign_change || direction_change
                    offset_idx = CURRENT_FRAME+round(j/oversampling_factor);
                    break;
                end
                offset_idx = CURRENT_FRAME+round(j/oversampling_factor);
            end
            
            % save saccade amplitude, saccade length, and seed location
            tmp_saccade_amp = pupil_position(offset_idx)-pupil_position(onset_idx);
%             if tmp_saccade_amp ~= 0
            if abs(tmp_saccade_amp) >= saccade_amplitude_threshold
                saccade(onset_idx,1) = tmp_saccade_amp;
                saccade_length(onset_idx,1) = offset_idx-onset_idx;
                saccade_seed(onset_idx,1) = CURRENT_FRAME - onset_idx; % for debugging. Indicates where the seed location (i.g., peak velocity) for each saccade relative to onset
            end

        end
    end
end

end
