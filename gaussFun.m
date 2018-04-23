function z = gaussFun(a,x0,y0,sx,sy,x,y,threshold)
    % Calculate value of a 2D gaussian
   
    if nargin < 8
        threshold = 0;
    end
    
    z = a.*exp(-((x-x0).^2./2./sx^2+(y-y0).^2./2./sy^2));
    
    % Use following logic to flatten top of Gaussian as desired
    if threshold
        ind = z>threshold;
        z(ind) = threshold;
    end
end