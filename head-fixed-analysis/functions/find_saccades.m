function [ saccade, saccade_length, velocity, saccade_seed ] = find_saccades(pupil_position, saccade_amplitude_threshold, saccade_velocity_threshold, saccade_onset_velocity_threshold,saccade_offset_velocity_threshold, frame_rate, plot_flag)
%find_saccades uses a velocity threshold to identify saccades. It ensures
%that a saccade has not been detected within a refractory period. It
%identifies saccade onset and offset as the first points in both
%directions from threshold crossing that go below a second threshold 
%(this must occur within a period of time defined by refractory_period). 

% INPUTS
% 
% pupil_position: Mx1 array of pupil position (deg)
% saccade_amplitude_threshold: (deg)
% saccade_velocity_threshold: (deg/sec)
% saccade_onset_velocity_threshold:  (deg/sec)
% saccade_offset_velocity_threshold:  (deg/sec)
% frame_rate: frames per second
% plot_flag: 0 (no plot) or 1 (plot trace with detected saccades)
%
% OUTPUTS
%
% saccade: Mx1 sparse array of saccade amplitudes. Values located at saccade onset frames.
% saccade_length: Mx1 sparse array of saccade lengths. Values located at saccade onset frames.
% velocity: Mx1 array of pupil velocity (deg/sec)
% saccade_seed: Mx1 sparse array of seed location (i.g., peak velocity) relative to saccade onset. Values located at saccade onset frames.


dt = 1/frame_rate; %time between frames in seconds
refractory_period = round(0.1/dt); % refractory period of 0.1s in frames

velocity = zeros(size(pupil_position));

%Calculate pupil velocity. This method eliminates the "jumps" from faulty
%tracking that incorrectly get detected as saccades. May want to address 
%this more concretely in the future 
for i = 2:length(pupil_position)-1
    velocity(i) = (pupil_position(i+1)-pupil_position(i-1))/(2*dt); % in seconds
end

velocityBool = abs(velocity)>saccade_velocity_threshold; %identify frames in which velocity exceeds threshold

saccade = zeros(length(velocity),1); %initialize saccade array
saccade_length = zeros(length(velocity),1); %initialize saccade length array
saccade_seed = zeros(length(velocity),1);

for i = refractory_period*2+1:length(velocityBool)-refractory_period-1
    if velocityBool(i) == 1 
        for j = -1:-1:-refractory_period
            if abs(velocity(i+j))<saccade_onset_velocity_threshold
                onset_idx = i+j;
                break;
            end
            onset_idx = i+j;
        end
        for j = 1:refractory_period
            if abs(velocity(i+j))<saccade_offset_velocity_threshold
                offset_idx = i+j;
                break;
            end
            offset_idx = i+j;
        end
        sacAmp = pupil_position(offset_idx)-pupil_position(onset_idx);
        if abs(sacAmp)>saccade_amplitude_threshold && ~any(saccade(onset_idx-refractory_period:onset_idx-1))
            saccade(onset_idx,1) = sacAmp;
            saccade_length(onset_idx,1) = offset_idx-onset_idx;
            saccade_seed(onset_idx,1) = i - onset_idx;
        end
    end
end

if plot_flag == 1
    figure;
    plot(pupil_position,'k')
    for i=1:length(saccade)
        if saccade(i)>0
            hold on;
            plot(i:i+saccade_length(i), pupil_position(i:i+saccade_length(i)),'r');
        end
        if saccade(i)<0
            hold on;
            plot(i:i+saccade_length(i), pupil_position(i:i+saccade_length(i)),'b');
        end
    end
    xlabel('Frame (50 fps)');
    ylabel('Azimuth (deg)');
end
end
