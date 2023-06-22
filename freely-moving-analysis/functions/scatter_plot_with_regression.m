function ax = scatter_plot_with_regression(x, y, varargin)
    
    % Parse inputs
    p = inputParser;
    addRequired(p, 'x');
    addRequired(p, 'y');
    addParameter(p, 'Size', 15);
    addParameter(p, 'Color', [0 0 0]);
    addParameter(p, 'XLim', [-100, 100]);
    addParameter(p, 'YLim', [-100, 100]);
    addParameter(p, 'PlotAxes', 'hv');
    addParameter(p, 'ax', false);
    parse(p, x, y, varargin{:})
    
    x = p.Results.x;
    y = p.Results.y;
    sz = p.Results.Size;
    c = p.Results.Color;
    xlimits = p.Results.XLim;
    ylimits = p.Results.YLim;
    PlotAxes = p.Results.PlotAxes;
    ax = p.Results.ax;
    
    % Plot data
    if ax == false
        figure; hold on; set(gcf, 'Position',  [200, 100, 200, 200]); hold on
        scatter(x, y, sz, c)
        ax = gca;
    else
        scatter(ax, x, y, sz, c)
    end

    xlim(xlimits)
    ylim(ylimits)

    if contains(PlotAxes, 'v')
        vline(0, 'k--')
    end
    if contains(PlotAxes, 'h')
        hline(0, 'k--')
    end

    % Regression line
    mdl = fitlm(x, y);
    plot(xlimits, predict(mdl, xlimits'), 'k-')
    
    % Regression stats
    text_x_offset = diff(xlimits)*0.05;
    text_y_offset = diff(ylimits)*0.05;
    text_spacing = diff(ylimits)*0.08;
    text(xlimits(1)+text_x_offset, ylimits(2) - text_y_offset - text_spacing*0, sprintf('R2 = %0.2f', mdl.Rsquared.Ordinary))
    text(xlimits(1)+text_x_offset, ylimits(2) - text_y_offset  - text_spacing*1, sprintf('Slope = %0.2f', mdl.Coefficients{2,1}))
    text(xlimits(1)+text_x_offset, ylimits(2) - text_y_offset  - text_spacing*2, sprintf('p = %0.3f', mdl.Coefficients{2,4}))

end