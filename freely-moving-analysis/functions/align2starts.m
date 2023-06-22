function [ startsMtrx, starts_trim ] = align2starts(sig, starts, window_size)
%ALIGNDATA Takes a continous time series as input and returns a matrix of
%arrays beggining at the indices listed in starts and 
%continuing for window_size time points. The dimensions of the returned 
%matrix will be (length(starts),window_size,number of input signals). If
%there are starts that do not fall within the signal, they will be
%discarded.
numSig = size(sig,2);
starts_trim = starts(starts>0 & (starts+window_size)<length(sig));
startsMtrx = zeros(length(starts_trim),window_size,numSig,'double');

for i = 1:numSig
    for j = 1:numel(starts_trim)
        startsMtrx(j,:,i) = sig(starts_trim(j):starts_trim(j)+window_size-1,i);
    end
    
%     hankelMtrx = hankel(sig(1:(end-window_size+1),i), sig((end-window_size+1):end,i)); % hankel is slow
%     startsMtrx(:,:,i) = hankelMtrx(starts_trim, :);
end
end

