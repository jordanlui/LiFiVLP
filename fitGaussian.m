function [model,x,y,z] = fitGaussian(x,y,z,k)
    % Fit a 2D Gaussian
    % Gaussian Function
    gaussFun = @(a,x0,y0,sx,sy,x,y) a .* exp(-( (x-x0).^2./2./sx^2 + (y-y0).^2./2./sy^2 ));
    x0 = [0.25 100 100 50 50]; % Initial Guess. [a x0 y0 sx sy]
    % Interpolate data, return meshgrid and vectors. 
    [xi,yi,zi,x,y,z] = prepGaussian(x,y,z,k);
    % Refine starting point estimate by looking at max z value
    [amp, ind] = max(z);
    x0 = [1 x(ind) y(ind) 50 50]; % Guess. [a x0 y0 sx sy]
    
    % Fit the model
    model = fit([x y],z,gaussFun,'Start',x0);
    
end