function Y = row_indexing(X,I)
% For each row in X return the value specified by the linear index I(row)
% This function pairs best with the function "first_event_after":
% Use first_event_after to find the indexes of the first event for each
% row, and then extract those events using row_indexing
%
% INPUTS
%
% X: MxN matrix (e.g., a matrix of trial data where each row is a trial)
% I: Mx1 array of indexes (e.g., frames at which you want extract data for each trial)
%
% OUTPUTS
%
% Y: Mx1 array values of X specified by I

if size(X,1) ~= numel(I)
    error('The number of rows in X must the number of values in I')
end

if size(I,1) < size(I,2)
    I = I';
end

if ismatrix(I)
    Y = NaN(size(I));
    for i = 1:numel(I)
        if ~isnan(I(i))
            Y(i) = X(i, I(i));
        else
            Y(i) = NaN;
        end
    end
end

end

