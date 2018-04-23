function ax = plotSparseSurface(x,y,z,dim,ax)
    % Plot sparse points as surface
    if nargin < 4
        dim = 100;
    end
    
    xmin=min(x);
    xmax=max(x);
    ymin=min(y);
    ymax=max(y);
    
    xi = linspace(xmin,xmax,dim);
    yi = linspace(ymin,ymax,dim);
    [xi,yi] = meshgrid(xi,yi);
    zi = griddata(x,y,z,xi,yi);
    
    ax = surf(xi,yi,zi);
%     contour(xi,yi,zi);
    xlabel('x'),ylabel('y'),zlabel('z')

end