function [bin_centers,counts] = violin_scatter(X,Y,c,scatter_size)
%violin_scatter Generate scatter plot of saccade endpoint (X) against initial position (Y) with histogram
%of endpoint (X) above. Colorcode scatter points by endpoint (X) coordinate.

if nargin == 3
    scatter_size = 36; % matlab default
elseif nargin == 4
    scatter_size = scatter_size;
end

% Generate saccade endpoint histogram 
bin_edges = linspace(-30,30,20); 
bin_centers = zeros(length(bin_edges)-1,1);
for i = 1:(length(bin_edges)-1)
    bin_centers(i) = mean(bin_edges(i:i+1));
end

counts = histcounts(X,linspace(-30,30,20),'Normalization','probability');
hold on; plot(bin_centers,counts,'color',c(size(c,1)/2,:)); hold on;

% set limits
ylim([-.3 0.5]);

% make some space under the density plot for the boxplot and raindrops
yl = get(gca, 'YLim');
set(gca, 'YLim', [-yl(2) yl(2)]);

% width of boxplot
wdth = yl(2) * 0.25;

% Determine drop Y position based on initial pupil position
drops_pos = Y/max(abs(Y))/5-yl(2)/2;
% drops_pos = normalize(Y,'scale',max(abs(Y)))/5-yl(2)/2;

% Determine color bin number based on initial pupil position
[bins, ~, ~] = createBins(Y,[min(Y) max(Y)],size(c,1));

% Create scatter
h = scatter(X,drops_pos, scatter_size, 'filled','MarkerEdgeColor','k', 'LineWidth', 0.1);

% Set color based on bin number
CData = zeros(size(Y,1),3);
for b = 1:size(c,1)
    CData_idx = find(bins==b);
    for i = 1:length(CData_idx)
        CData(CData_idx(i),:) = c(b,:);
    end
end
h.CData = CData;
end

