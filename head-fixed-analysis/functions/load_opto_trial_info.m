function [tmp_opto_duration,tmp_opto_intensity,tmp_opto_cyclelength,tmp_opto_pulselength] = load_opto_trial_info(puffData,experiment)
%load_opto_trial_info Combines opto onset data from puffData struct (nidaq data)
%with opto intensity,duration,cyclelength,and pulselength info from experiment
%struct (matlab data)
opto_trial_idx = find(experiment.trials.opto==1); 
num_completed_opto_trials = sum(puffData.optoTrial);
completed_opto_trial_idx = opto_trial_idx(1:num_completed_opto_trials);

tmp_opto_duration = puffData.optoTrial;
tmp_opto_intensity = puffData.optoTrial;
tmp_opto_cyclelength = puffData.optoTrial;
tmp_opto_pulselength = puffData.optoTrial;

tmp_opto_duration(tmp_opto_duration==1) = experiment.trials.opto_duration(completed_opto_trial_idx);
tmp_opto_intensity(tmp_opto_intensity==1) = experiment.trials.opto_intensity(completed_opto_trial_idx);
tmp_opto_cyclelength(tmp_opto_cyclelength==1) = experiment.trials.opto_cyclelength(completed_opto_trial_idx);
tmp_opto_pulselength(tmp_opto_pulselength==1) = experiment.trials.opto_pulselength(completed_opto_trial_idx);
end

