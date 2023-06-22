function [new_response_data, new_training_data] = downsample_majority_class(response, predictors)
% downsample_majority_class resamples response and predictors so that there are an equal number
% of values for each category

% get counds of values for each category
[C,~,ic] = unique(response);
category_counts = accumarray(ic,1);

% use minimum count for downsampling
numToSample = min(category_counts);

% randomly resample each class so that so that there are equal numbers of values
new_response_data = [];
new_training_data = [];
for i = 1:numel(C)
    tmp_sample_idx = randsample(category_counts(i), numToSample);
    
    tmp_category_response_data = response(response == C(i));
    new_response_data = [new_response_data; tmp_category_response_data(tmp_sample_idx)];
    
    tmp_category_training_data = predictors(response == C(i), :);
    new_training_data = [new_training_data; tmp_category_training_data(tmp_sample_idx, :)];
    
end

end

