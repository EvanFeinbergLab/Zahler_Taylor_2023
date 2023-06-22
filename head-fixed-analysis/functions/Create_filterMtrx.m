function [filterMtrx] = Create_filterMtrx(hankelMtrx,FILTER)
%Create_filterMtrx takes a hankel matrix constructed from all starts
%throughout the combined sessions and filters it using the boolean FILTER

fprintf('\nFiltering %d of %d trials\n', sum(FILTER), numel(FILTER));

% Filters all fields in hankelMtrx
fields = fieldnames(hankelMtrx);
for i = 1:numel(fields)
    filterMtrx.(fields{i}) = hankelMtrx.(fields{i})(FILTER,:);
end

end