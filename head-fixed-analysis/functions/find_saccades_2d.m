function [ saccadeAmp, saccadeAngle, saccadeX, saccadeY, saccade_length ] = find_saccades_2d(pupilX, pupilY, saccade_amplitude_threshold, saccade_velocity_threshold, saccade_onset_velocity_threshold, saccade_offset_velocity_threshold, frame_rate, plot_flag)
%FIND_SACCADES_2D
% 1) Velocity threshold ("saccade_velocity_threshold") is used to identify candidate saccades
% 2) Saccade candidates are discarded if they fall within a 100ms refractory period after other candidate saccades
% 3) Saccade onset is defined as the point before the velocity threshold
% crossing where the velocity dips below a second threshold ("saccade_onset_velocity_threshold")
% 4) Saccade offset is the point after the velocity threshold crossing where velocity dips below a
% third threshold ("saccade_offset_velocity_threshold").
% 5) Finally, saccade amplitude must exceed a final threshold ("saccade_amplitude_threshold")

% Velocity thresholds should be in deg/sec
% Amplitude threshold should be in deg

dt = 1/frame_rate; %time between frames in seconds
refractory_period = round(0.1/dt); % refractory period of 0.1s in frames

%calculate 2d velocity (in seconds)
velocity = zeros(size(pupilX));
for i = 2:length(pupilX)-1
    velocity(i) = sqrt((pupilX(i+1)-pupilX(i-1))^2 + (pupilY(i+1)-pupilY(i-1))^2)/(2*dt);
end

%identify frames in which velocity exceeds threshold
velocityBool = abs(velocity)>saccade_velocity_threshold;

% initialize output arrays
saccadeAmp = zeros(length(velocity),1); % initialize saccade 2d amplitude array
saccadeAngle = NaN(length(velocity),1); % initialize saccade angle array
saccadeX = zeros(length(velocity),1); % initialize saccade x component array
saccadeY = zeros(length(velocity),1); % initialize saccade y component array
saccade_length = NaN(length(velocity),1); % initialize saccade duration array

% 
for i = refractory_period+1:length(velocityBool)-refractory_period-1
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

        [theta, rho] = cart2pol(pupilX(offset_idx)-pupilX(onset_idx), pupilY(offset_idx)-pupilY(onset_idx));
        sacAmp2D = rho;
        sacAmpX = pupilX(offset_idx)-pupilX(onset_idx);
        sacAmpY = pupilY(offset_idx)-pupilY(onset_idx);
        
        if onset_idx > refractory_period
            if abs(sacAmp2D)>saccade_amplitude_threshold  && ~any(saccadeAmp(onset_idx-refractory_period:onset_idx-1))
                saccadeAmp(onset_idx,1) = rho;
                saccadeAngle(onset_idx,1) = mod(rad2deg(theta), 360);
                saccadeX(onset_idx,1) = sacAmpX;
                saccadeY(onset_idx,1) = sacAmpY;
                saccade_length(onset_idx,1) = offset_idx-onset_idx;
            end
        end
        
    end
end

if plot_flag == 1
    figure;
    ax1 = subplot(2,1,1); hold on; title('Pupil X')
    plot(pupilX,'k')
    for i=1:length(saccadeAmp)
        if saccadeAmp(i)>0
            hold on;
            plot(i:i+saccade_length(i), pupilX(i:i+saccade_length(i)),'r', 'LineWidth', 3);
        end
        if saccadeAmp(i)<0
            hold on;
            plot(i:i+saccade_length(i), pupilX(i:i+saccade_length(i)),'b', 'LineWidth', 3);
        end
    end
    xlabel('Frame (50 fps)');
    ylabel('Azimuth (deg)');
    
    
    ax2 = subplot(2,1,2); hold on; title('Pupil Y')
    plot(pupilY,'k')
    for i=1:length(saccadeAmp)
        if saccadeAmp(i)>0
            hold on;
            plot(i:i+saccade_length(i), pupilY(i:i+saccade_length(i)),'r', 'LineWidth', 3);
        end
        if saccadeAmp(i)<0
            hold on;
            plot(i:i+saccade_length(i), pupilY(i:i+saccade_length(i)),'b', 'LineWidth', 3);
        end
    end
    xlabel('Frame (50 fps)');
    ylabel('Altitude (deg)');
    
    linkaxes([ax1, ax2])
    
%     figure;
%     theta = deg2rad(0:15:360);
%     rho = histcounts(saccadeAngle, 0:15:360);
%     rho = [rho rho(1)];
%     polarplot(theta, rho)

    
end
end
