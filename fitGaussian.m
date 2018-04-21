function [model,x,y,z] = fitGaussian(x,y,z,k,symmetric)
    % Fit a 2D Gaussian
    % Gaussian Function
    posTol = 800; % Tolerance for spatial value
    sigma0 = 300; % Guess for Gaussian width
    sigmaTol = 1000;
    
    if nargin < 5
        % If not specified, we will fit an asymmetric Gaussian
        symmetric = 0;
    end
    
    % Interpolate data, return meshgrid and vectors. 
    [xi,yi,zi,x,y,z] = prepGaussian(x,y,z,k);
    % Refine starting point estimate by looking at max z value
    [amp, ind] = max(z);
    
    if symmetric
        gaussFun = @(a,x0,y0,sx,x,y) a .* exp(-( (x-x0).^2./2./sx^2 + (y-y0).^2./2./sx^2 ));
        x0 = [amp x(ind) y(ind) sigma0]; % Guess. [a x0 y0 sx sy]
        xU = [5 posTol posTol sigmaTol];
        xL = [0 -posTol -posTol 0];
    else
        gaussFun = @(a,x0,y0,sx,sy,x,y) a .* exp(-( (x-x0).^2./2./sx^2 + (y-y0).^2./2./sy^2 ));
        x0 = [amp x(ind) y(ind) sigma0 sigma0]; % Guess. [a x0 y0 sx sy]
        xU = [5 posTol posTol sigmaTol sigmaTol];
        xL = [0 -posTol -posTol 0 0];
    end
    
    % Fit the model
    model = fit([x y],z,gaussFun,'Start',x0,...
        'Lower',xL,...
        'Upper',xU);
    
end