function [trainedClassifier] = train_multinomial_classifier(training_data, response_data)
%  Input:
%      training_data: 
%
%      response_data: ... The length of response_data and the number of
%       rows of training_data must be equal.
%
%  Output:
%      trained_classifier: A struct containing the trained classifier.
%
%      trained_classifier.predict_fcn: A function to make predictions on new
%       data.
%
%      validation_accuracy: A double containing the accuracy in percent.
%
% Use the code to train the model with new data. To retrain your
% classifier, call the function from the command line with your original
% data or new data as the input arguments trainingData and responseData.
%
% For example, to retrain a classifier trained with the original data set T
% and response Y, enter:
%   [trained_classifier, validation_accuracy] = train_logistic_classifier(T, Y)
%
% To make predictions with the returned 'trainedClassifier' on new data T2,
% use
%   yfit = trainedClassifier.predictFcn(T2)
%

% Train using fitglm
[B,dev,stats] = mnrfit(...
    training_data, ...
    response_data);

% [B,dev,stats] = mnrfit(...
%     training_data, ...
%     response_data, ...
%     'model','ordinal');


% Convert predicted probabilities to predicted class labels and scores.

findIdxOfMax = @(p) find(max(p) == p);
convertSuccessProbsToPredictions =  @(p) cellfun(findIdxOfMax, num2cell(p, 2));
returnMultipleValuesFcn = @(varargin) varargin{1:max(1,nargout)};
predictionsAndProbsFcn = @(p) returnMultipleValuesFcn( convertSuccessProbsToPredictions(p), p );

removeNaN = @(p) fillmissing(p, 'constant', 1);

% Create the result struct with predict function
multinomialLogisticRegressionPredictFcn = @(x) predictionsAndProbsFcn( removeNaN(mnrval(B, x)) );
trainedClassifier.predictFcn = @(x) multinomialLogisticRegressionPredictFcn(x);

% Add additional fields to the result struct
trainedClassifier.model.B = B;
trainedClassifier.model.dev = dev;
trainedClassifier.model.stats = stats;
