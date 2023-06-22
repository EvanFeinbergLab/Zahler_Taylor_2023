function [ mtrx ] = combine_sessions(BehaviorFiles, trial_param, normalize_sessions, varargin)
%combine_sessions Combine data across sessions

% Optional name-value inputs:
% 'FindSaccadesEye': 'mean' (default), 'left', or 'right' determines which eye trace input to find_saccades and find_saccades_2d functions

% Parse other inputs
p = inputParser;

addRequired(p,'BehaviorFiles');
addRequired(p,'trial_param');
addRequired(p,'normalize_sessions');
addParameter(p,'FindSaccadesEye', 'mean');
addParameter(p,'ZScoreStrainGauge', true);

parse(p, BehaviorFiles, trial_param, normalize_sessions, varargin{:})

BehaviorFiles = p.Results.BehaviorFiles;
trial_param = p.Results.trial_param;
normalize_sessions = p.Results.normalize_sessions;
FindSaccadesEye = p.Results.FindSaccadesEye;
ZScoreStrainGauge = p.Results.ZScoreStrainGauge;

%Initialize variables that you are combining
leftPuff=[];
rightPuff = [];
leftPuff2=[];
rightPuff2 = [];
pupilX_L = [];
pupilY_L=[];
pupilX_R=[];
pupilY_R=[];
pupilX_Mean = [];
pupilY_Mean = [];
eyelid_L=[];
eyelid_R=[];
pupilArea_L=[];
pupilArea_R=[];
strain_gauge = []; 
optoTrial = [];
opto_duration = [];
opto_intensity = [];
opto_cyclelength = [];
opto_pulselength = [];
mouseID = [];

BehaviorFiles_mouse_id = {};
for i = 1:length(BehaviorFiles)
    BehaviorFiles_mouse_id{i} = BehaviorFiles(i).config.experiment.mouse_id;
end
unique_mouse_ids = unique(BehaviorFiles_mouse_id);


for i = 1:length(BehaviorFiles)
    
    %load session variables 
    experiment = BehaviorFiles(i).experiment;
    puffData = BehaviorFiles(i).puffData;
    
    % correct some sessions that have an erroneous timestamp at first frame
    puffData.leftPuff(1) = 0;
    puffData.rightPuff(1) = 0;
    puffData.optoTrial(1) = 0;
    if isfield(puffData, 'leftPuff2')
        puffData.leftPuff2(1) = 0;
        puffData.rightPuff2(1) = 0;
    end
    
    
    dlcPupilData = BehaviorFiles(i).dlcPupilData;
    dlcLikelihoodData = BehaviorFiles(i).dlcLikelihoodData;
    dlcCamera2Data = BehaviorFiles(i).dlcCamera2Data;
    dlcCamera2LikelihoodData = BehaviorFiles(i).dlcCamera2LikelihoodData;
    tmp_strain_gauge = BehaviorFiles(i).strain_gauge;
    % If strain_gauge velocity exceeds 2 (200 mV/s), replace with NaN 
    % b/c it is an artifact of the DAC reaching its threshold and resetting
    tmp_strain_gauge(abs(diff(tmp_strain_gauge))>2) = NaN; 
    
    % encode mouse_id as a number
    
    tmp_mouseID = ones(size(dlcPupilData,1),1)*find(strcmp(BehaviorFiles(i).config.experiment.mouse_id, unique_mouse_ids));
    mouseID = vertcat(mouseID, tmp_mouseID); 

    %load and combine opto trial info across sessions
    [tmp_opto_duration,tmp_opto_intensity,tmp_opto_cyclelength,tmp_opto_pulselength] = load_opto_trial_info(puffData,experiment);
    optoTrial = vertcat(optoTrial, puffData.optoTrial); 
    opto_duration = vertcat(opto_duration, tmp_opto_duration);
    opto_intensity = vertcat(opto_intensity, tmp_opto_intensity);
    opto_cyclelength = vertcat(opto_cyclelength, tmp_opto_cyclelength);
    opto_pulselength = vertcat(opto_pulselength, tmp_opto_pulselength);            
    
    % Compute pupil position in degrees using Sakatani and Isa method. 
    % Prior to 2021-09-27, lens mag was specific to each setup
    if str2double(BehaviorFiles(i).config.experiment.date) < 20210927
        if contains(BehaviorFiles(i).config.computer_name, 'Sebi 2P Room')
            lens_mag = 104; % lens mag = pixels/mm
        else
            lens_mag = 62; % Sebi's old lenses
        end
    else
        lens_mag = 104; % lens mag = pixels/mm
    end
    
    %Left eye 
    [tmpPupilX_L_deg, tmpPupilY_L_deg, tmpEyelid_L, tmpPupilArea_L] = extract_pupil_trace(dlcPupilData, dlcLikelihoodData, lens_mag, 'sakatani');
    
    %Right eye
    [tmpPupilX_R_deg, tmpPupilY_R_deg, tmpEyelid_R, tmpPupilArea_R] = extract_pupil_trace(dlcCamera2Data, dlcCamera2LikelihoodData, lens_mag, 'sakatani');

    if normalize_sessions == 0 
        pupilX_L = vertcat(pupilX_L, tmpPupilX_L_deg);
        pupilY_L = vertcat(pupilY_L, tmpPupilY_L_deg);
        pupilX_R = vertcat(pupilX_R, tmpPupilX_R_deg);
        pupilY_R = vertcat(pupilY_R, tmpPupilY_R_deg);
        strain_gauge = vertcat(strain_gauge, tmp_strain_gauge-nanmean(tmp_strain_gauge)); %subtract mean for individual sessions because of differences in strain gauge zeroing
    elseif normalize_sessions == 1
        pupilX_L = vertcat(pupilX_L, tmpPupilX_L_deg-nanmean(tmpPupilX_L_deg));
        pupilY_L = vertcat(pupilY_L, tmpPupilY_L_deg-nanmean(tmpPupilY_L_deg));
        pupilX_R = vertcat(pupilX_R, tmpPupilX_R_deg-nanmean(tmpPupilX_R_deg));
        pupilY_R = vertcat(pupilY_R, tmpPupilY_R_deg-nanmean(tmpPupilY_R_deg));
        pupilX_Mean = vertcat(pupilX_Mean, (tmpPupilX_L_deg+tmpPupilX_R_deg)./2-nanmean((tmpPupilX_L_deg+tmpPupilX_R_deg)./2));
        pupilY_Mean = vertcat(pupilY_Mean, (tmpPupilY_L_deg+tmpPupilY_R_deg)./2-nanmean((tmpPupilY_L_deg+tmpPupilY_R_deg)./2));
        
        %Z-Score strain_gauge 
        if ZScoreStrainGauge == true
            strain_gauge = vertcat(strain_gauge, (tmp_strain_gauge-nanmean(tmp_strain_gauge))./nanstd(tmp_strain_gauge)); 
        elseif ZScoreStrainGauge == false
            strain_gauge = vertcat(strain_gauge, tmp_strain_gauge-nanmean(tmp_strain_gauge)); 
        end
    end
    eyelid_L = vertcat(eyelid_L,(tmpEyelid_L-nanmean(tmpEyelid_L))./nanstd(tmpEyelid_L));
    eyelid_R = vertcat(eyelid_R,(tmpEyelid_R-nanmean(tmpEyelid_R))./nanstd(tmpEyelid_R));
    pupilArea_L = vertcat(pupilArea_L,(tmpPupilArea_L-nanmean(tmpPupilArea_L))./nanstd(tmpPupilArea_L)); %zscore before concatenating because of differing light conditions across days
    pupilArea_R = vertcat(pupilArea_R,(tmpPupilArea_R-nanmean(tmpPupilArea_R))./nanstd(tmpPupilArea_R));
    leftPuff = vertcat(leftPuff, puffData.leftPuff);
    rightPuff = vertcat(rightPuff, puffData.rightPuff);
    
    if isfield(puffData, 'leftPuff2')
        leftPuff2 = vertcat(leftPuff2, puffData.leftPuff2);
        rightPuff2 = vertcat(rightPuff2, puffData.rightPuff2);
    else
        leftPuff2 = vertcat(leftPuff2, zeros(size(puffData.leftPuff)));
        rightPuff2 = vertcat(rightPuff2, zeros(size(puffData.rightPuff)));
    end
   
    
end

if normalize_sessions == 0
    % subtract mean from left and right pupil positions and store
    % separately
    mtrx.pupilX_L = pupilX_L-nanmean(pupilX_L); 
    mtrx.pupilY_L = pupilY_L-nanmean(pupilY_L);
    mtrx.pupilX_R = pupilX_R-nanmean(pupilX_R);
    mtrx.pupilY_R = pupilY_R-nanmean(pupilY_R);
    % Average left and right pupil positions and subtract mean
    mtrx.pupilX_Mean = (pupilX_L+pupilX_R)./2-nanmean((pupilX_L+pupilX_R)./2);
    mtrx.pupilY_Mean = (pupilY_L+pupilY_R)./2-nanmean((pupilY_L+pupilY_R)./2);
    %Z-Score strain_gauge 
    if ZScoreStrainGauge == true
        mtrx.strain_gauge = strain_gauge./nanstd(strain_gauge);
    end
elseif normalize_sessions == 1
    %Individual sessions already normalized
    mtrx.pupilX_L = pupilX_L;
    mtrx.pupilY_L = pupilY_L;
    mtrx.pupilX_R = pupilX_R;
    mtrx.pupilY_R = pupilY_R;
    mtrx.pupilX_Mean = pupilX_Mean;
    mtrx.pupilY_Mean = pupilY_Mean;
    mtrx.strain_gauge = strain_gauge;
end
mtrx.pupilArea = (pupilArea_L+pupilArea_R)./2;
mtrx.eyelid = (eyelid_L+eyelid_R)./2;
mtrx.leftPuff = leftPuff;
mtrx.rightPuff = rightPuff;
mtrx.leftPuff2 = leftPuff2;
mtrx.rightPuff2 = rightPuff2;
mtrx.allPuff = leftPuff+rightPuff+leftPuff2+rightPuff2;

saccade_onset_velocity_threshold = 30; % degrees/second
saccade_offset_velocity_threshold = 20; % degrees/second
switch FindSaccadesEye
    case 'mean'
        [mtrx.saccade, mtrx.saccade_length] = find_saccades(mtrx.pupilX_Mean,trial_param.saccade_threshold,trial_param.velocity_threshold,saccade_onset_velocity_threshold,saccade_offset_velocity_threshold,trial_param.frame_rate,0);
        [mtrx.saccade_2d_amp, mtrx.saccade_2d_angle, mtrx.saccade_2d_x, mtrx.saccade_2d_y ] = find_saccades_2d(mtrx.pupilX_Mean, mtrx.pupilY_Mean, trial_param.saccade_threshold,trial_param.velocity_threshold,saccade_onset_velocity_threshold,saccade_offset_velocity_threshold,trial_param.frame_rate,0);
    case 'left'
        [mtrx.saccade, mtrx.saccade_length] = find_saccades(mtrx.pupilX_L,trial_param.saccade_threshold,trial_param.velocity_threshold,saccade_onset_velocity_threshold,saccade_offset_velocity_threshold,trial_param.frame_rate,0);
        [mtrx.saccade_2d_amp, mtrx.saccade_2d_angle, mtrx.saccade_2d_x, mtrx.saccade_2d_y ] = find_saccades_2d(mtrx.pupilX_L, mtrx.pupilY_L, trial_param.saccade_threshold,trial_param.velocity_threshold,saccade_onset_velocity_threshold,saccade_offset_velocity_threshold,trial_param.frame_rate,0);
    case 'right'
        [mtrx.saccade, mtrx.saccade_length] = find_saccades(mtrx.pupilX_R,trial_param.saccade_threshold,trial_param.velocity_threshold,saccade_onset_velocity_threshold,saccade_offset_velocity_threshold,trial_param.frame_rate,0);
        [mtrx.saccade_2d_amp, mtrx.saccade_2d_angle, mtrx.saccade_2d_x, mtrx.saccade_2d_y ] = find_saccades_2d(mtrx.pupilX_R, mtrx.pupilY_R, trial_param.saccade_threshold,trial_param.velocity_threshold,saccade_onset_velocity_threshold,saccade_offset_velocity_threshold,trial_param.frame_rate,0);
end
mtrx.optoTrial = optoTrial;
mtrx.opto_duration = opto_duration;
mtrx.opto_intensity = opto_intensity;
mtrx.opto_cyclelength = opto_cyclelength;
mtrx.opto_pulselength = opto_pulselength;
mtrx.mouseID = mouseID;
end

