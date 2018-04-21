function plotCDF(y)
    % Plot CDF
    
    % Option 1: Fit a Distribution
%     x = linspace(min(y),max(y));
%     pd = fitdist(y,'Normal');
%     y = cdf(pd,x);
%     plot(x,y);
    
    % Option 2: Discrete CDF, don't fit disn
    x = sort(y);
    y = cumsum(y) ./ max(y);
    plot(x,y);
    ylim([-10, 110])
    xlim([-10, 10+round(length(y),-1)])
    
    % Label result
    xlabel('error (mm)'), ylabel('Fraction of total')
    title('Cumulative Distribution Fucntion')
    
end