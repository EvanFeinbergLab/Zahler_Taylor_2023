function [first_event_idx, first_event_logical] = first_event_after(event_matrix, after)
% Find the first non-zero value after a given timepoint. Useful for
% determining at what timepoint a saccade or head movement occured, but can
% be used for other situations.
%
% INPUTS
% event_matrix: MxN matrix representing trials (each row as one trial)
% after: timepoint (frame) after which you want to find a non-zero value in event_matrix (includes the value represented by "after")
%
% OUTPUTS
% first_event_idx: Mx1 array with the index (starting from the trial start)
% of the first event for each trial. If there was no event, the entry will
% be NaN
%
% first_event_logical: same information as first_event_idx, but data is
% represented as an MxN matrix (same dimensions as event_matrix)

[~, I_logical] = find_2d(event_matrix, [], 2);
I_logical(:,1:after-1) = false;
[first_event_idx, first_event_logical] = find_2d(I_logical, 1, 2);
first_event_idx = cell2mat(first_event_idx);

end

