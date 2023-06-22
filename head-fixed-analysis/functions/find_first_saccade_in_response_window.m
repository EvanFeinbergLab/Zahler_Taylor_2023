function [response_window_first_saccade, response_window_first_saccade_index] = find_first_saccade_in_response_window(saccade,response_window)
%find_first_saccade_in_response_window Finds the first saccade in the 
%response window. Returns amplitudes and indices. Returns NaN if no saccade
%in response window
%   saccade: # trials x # timepoints size matrix with saccade amplitudes at
%   saccade onset and 0 otherwise
%   response_window: timepoints within which to identify the first saccade

    response_window_first_saccade = zeros(size(saccade,1),1);
    response_window_first_saccade_index = zeros(size(saccade,1),1);
    for i = 1:size(saccade,1)
        response_window_saccade_indices = find(saccade(i,response_window)~=0);
        if ~isempty(response_window_saccade_indices)
            response_window_first_saccade_index(i) = response_window_saccade_indices(1);
            response_window_first_saccade(i) = saccade(i,response_window(response_window_saccade_indices(1)));
        else
            response_window_first_saccade_index(i) = NaN;
            response_window_first_saccade(i) = NaN;
        end
    end
end

