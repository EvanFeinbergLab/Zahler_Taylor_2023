function [bins, bin_array, bin_centers] = createBins( signal,bin_edges, num_bins )
%createBins takes a signal and the number of bins you'd like to split it into. 
%Generates an array with bin limits such that each bin contains the same
%number of data points
%data points.

bin_array = linspace(bin_edges(1),bin_edges(2),num_bins+1);
bins = discretize(signal, bin_array);

bin_centers = zeros(length(bin_array)-1,1);
for i = 1:(length(bin_array)-1)
    bin_centers(i) = mean(bin_array(i:i+1));
end

end

