function [I, I_logical] = find_2d(X, K, dim)
% Similar to "find" function, except that it can operate on arrays.
%
% E.g., find_2d(X, K, 1) will apply "find" to each column
% Similarly, find_2d(X, K, 2) will apply "find" to each row
%
% I is a cell array of indices for each column or row
% I_logical is a matrix with the same dimensions as X where each non-zero
% value of X is replaced with a 1

X_size = size(X);

I_logical = logical(zeros(size(X)));

if dim == 1
    I = cell(1,X_size(2));
    for i = 1:X_size(2)
        if ~isempty(K)
            I{i} = find(X(:,i), K);
            if ~isempty(I{i})
                I_logical(I{i},i) = 1;
            end
        else
            I{i} = find(X(:,i));
            if ~isempty(I{i})
                I_logical(I{i},i) = 1;
            end
        end
    end

elseif dim == 2
    I = cell(X_size(1),1);
    for j = 1:X_size(1)
        if ~isempty(K)
            I{j} = find(X(j,:), K);
            if ~isempty(I{j})
                I_logical(j, I{j}) = 1;
            end
        else
            I{j} = find(X(j,:));
            if ~isempty(I{j})
                I_logical(j, I{j}) = 1;
            end
        end
    end
end

I(  cellfun(@isempty, I)  ) = {NaN};


end

