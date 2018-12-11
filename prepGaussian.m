function [xi,yi,zi,x,y,z] = prepGaussian(x,y,z,k)
    % Interpolate data to fit meshgrid and fit the Gaussian via Matlab
    x = x(1:k:end);
    y = y(1:k:end);
    z = z(1:k:end);
    
    xmin=min(x);
    xmax=max(x);
    ymin=min(y);
    ymax=max(y);
    dim = 51;
    
    xi = linspace(xmin,xmax,dim);
    yi = linspace(ymin,ymax,dim);
    [xi,yi] = meshgrid(xi,yi);
    zi = griddata(x,y,z,xi,yi);
    
    % Cut out NAN values and return the vectors
    x = reshape(xi,[numel(xi),1]);
    y = reshape(yi,[numel(yi),1]);
    z = reshape(zi,[numel(zi),1]);
    t = ~isnan(x) & ~isnan(y) & ~isnan(z);
    x=x(t); y=y(t); z=z(t);

end