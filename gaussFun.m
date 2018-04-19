function z = gaussFun(a,x0,y0,sx,sy,x,y)
    % Calculate value of a 2D gaussian
    
    z = a.*exp(-((x-x0).^2./2./sx^2+(y-y0).^2./2./sy^2));
end