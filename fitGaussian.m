function [model,x,y,z] = fitGaussian(x,y,z,k,symmetric,threshold)
    % Fit a 2D Gaussian
    % Note this code needs much improving. Currently cannot fit symmetric
    % thresholded Gaussian
    
    posTol = 800; % Tolerance for spatial value
    sigmaTol = 1000; % Tolerance for size of sigma value
    sigma0 = 300; % Guess for Gaussian width
    threshTol = 2; % Max tolerance for Gaussian cutoff height
    threshold0 = 0.3;
    
    if nargin < 5
        % If not specified, we will fit an asymmetric Gaussian with no
        % threshold
        symmetric = 0;
        threshold = 0;
    end
    
    % Interpolate data, return meshgrid and vectors. 
    [xi,yi,zi,x,y,z] = prepGaussian(x,y,z,k);
    % Refine starting point estimate by looking at max z value
    [amp, ind] = max(z);
    % Set tolerances
    x0 = [amp x(ind) y(ind) sigma0 sigma0]; % Guess. [a x0 y0 sx sy]
    xU = [5 posTol posTol sigmaTol sigmaTol];
    xL = [0 -posTol -posTol 0 0];
    
    if symmetric
        f = @(a,x0,y0,sx,x,y) a .* exp(-( (x-x0).^2./2./sx^2 + (y-y0).^2./2./sx^2 ));
        x0 = x0(1:4);
        xU = xU(1:4);
        xL = xL(1:4);
    else
        f = @(a,x0,y0,sx,sy,x,y) a .* exp(-( (x-x0).^2./2./sx^2 + (y-y0).^2./2./sy^2 ));
        
    end
    
    if threshold
        if threshold ~=1 % If specific guess was supplied, load it
            threshold0 = threshold;
        end
        f = @(a,x0,y0,sx,sy,threshold,x,y) gaussFun(a,x0,y0,sx,sy,x,y,threshold);
        x0 = [x0 threshold0];
        xU = [xU threshTol];
        xL = [xL 0];
    end
    
    % Fit the model
    model = fit([x y],z,f,'Start',x0,...
        'Lower',xL,...
        'Upper',xU);
    
end