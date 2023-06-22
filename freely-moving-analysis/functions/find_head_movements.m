function [head_movement_amplitude, head_movement_length, velocity, head_movement_seed] = find_head_movements(head_trace, frame_rate, head_movement_velocity_threshold, head_movement_minimum_duration)
%FIND_HEAD_MOVEMENTS Identifies the start, end, and amplitude of head movements
% Adapted from Masullo et al., 2019
%
% Head movement events are identified when angular head movement velocity exceed head_movement_velocity_threshold (°/s)
% in the same direction for at least head_movement_minimum_duration (frames). Head movement onset and offset are found by
% moving backward or forward, respectively, in time from the seed location until velocity fell below 20% of vel_threshold
% or velocity direction changes. Head movement amplitude is the difference between the head positions at offset and onset.
%
%
% INPUTS
%
% head_trace: Mx1 array of yaw, pitch, or roll head movement positions in degrees
% frame_rate: data framerate (frames per second)
% head_movement_velocity_threshold: (deg/s)
% head_movement_minimum_duration: (frames)
%
%
% OUTPUTS
%
% head_movement_amplitude: Mx1 sparse array of head movement amplitudes. Values located at saccade onset frames.
% head_movement_length: Mx1 sparse array of saccade lengths. Values located at saccade onset frames.
% velocity: Mx1 sparse array of seed location (i.g., peak velocity) relative to saccade onset. Values located at saccade onset frames.


% invert pupil_position if M<N
if size(head_trace,1) < size(head_trace,2) || size(head_trace,2) ~= 1
    head_trace = head_trace';
end

% calculate velocity
dt = 1/frame_rate; %time between frames in seconds
velocity = [0; diff(head_trace)]/dt;

%initialize head movement arrays
head_movement_amplitude = zeros(size(velocity));
head_movement_length = zeros(size(velocity));
head_movement_seed = zeros(length(velocity),1);

previous_movement_end = 0;
for CURRENT_FRAME = head_movement_minimum_duration:numel(head_trace)

    % Velocity must exceed head_movement_velocity_threshold and remain in the same direction for a period of time (head_movement_minimum_duration)
    if all(abs(velocity(CURRENT_FRAME-head_movement_minimum_duration+1:CURRENT_FRAME)) > head_movement_velocity_threshold) && numel(unique(sign(velocity(CURRENT_FRAME-head_movement_minimum_duration+1:CURRENT_FRAME))))==1 
        
        % find movement onset
        tmp_movement_start = CURRENT_FRAME - head_movement_minimum_duration+1;
        if tmp_movement_start > 1
            for j = 1:1:100
                % keep advancing backward in time until velocity sign switches or velocity falls to 20% of head_movement_velocity_threshold
                if tmp_movement_start - j > 1 && sign(velocity(tmp_movement_start - j)) == sign(velocity(tmp_movement_start)) && abs(velocity(tmp_movement_start - j)) > head_movement_velocity_threshold/5
                    continue;
                else
                    break;
                end
            end
            tmp_movement_start = tmp_movement_start - j;
        end
        
        % find movement offset
        tmp_movement_end = CURRENT_FRAME;
        if tmp_movement_end < numel(velocity)
            for j = 1:1:100
                % keep advancing forward in time until velocity sign switches or velocity falls to 20% of head_movement_velocity_threshold
                if tmp_movement_end + j < numel(velocity) && sign(velocity(tmp_movement_end + j)) == sign(velocity(tmp_movement_end)) && abs(velocity(tmp_movement_end + j)) > head_movement_velocity_threshold/5
                    continue;
                else
                    break;
                end
            end
            tmp_movement_end = tmp_movement_end + j;
        end
        
        % make sure there's no overlab between movements and save movement amplitude/length 
        if tmp_movement_start > previous_movement_end
            head_movement_amplitude(tmp_movement_start) = head_trace(tmp_movement_end) - head_trace(tmp_movement_start);
            head_movement_length(tmp_movement_start) = tmp_movement_end - tmp_movement_start;
            head_movement_seed(tmp_movement_start) = CURRENT_FRAME - tmp_movement_start; % for debugging. Indicates where the seed location for each head movement relative to onset
            previous_movement_end = tmp_movement_end;
        end
        
    end
    
end

end

