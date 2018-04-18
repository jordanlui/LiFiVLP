function [f1, fitresult] = compareGaussianFit(x,y,z,x0,gaussFun,grid)
    % Plot 1: Fit data with two different Gaussians and plot to compare
    figure()
    subplot(2,2,1)
    hold on
    plot3(x,y,z,'r.'),xlabel('x'),ylabel('y'), title('original data, no crop')
    contour(grid.xi,grid.yi,grid.zi)
    
      
    % Plot 2: Display cropped data
    t = x>1 & y>1 & x<350 & y<350;
    % t = ones(length(x),1); % Alternative if we want to fit all data
    subplot(2,2,2)
    plot3(x(t),y(t),z(t),'r.'),xlabel('x'),ylabel('y'), title('Cropped Data')

    % Plot 3: Gaussian fit with system default equations
    % Plot with no clipping
    f1 = fit([x y],z,gaussFun,'Start',x0);
    zpred = gaussFun(f1.a,f1.x0,f1.y0,f1.sx,f1.sy,x,y);
    error = rms(zpred-z);
    subplot(2,2,3)
    plot(f1,[x,y],z),xlabel('x'),ylabel('y'), title(sprintf('Gaussian fit. RMS:%.2f',error))
    
    % Plot with clipping
%     f1 = fit([x(t) y(t)],z(t),gaussFun,'Start',x0);
%     subplot(2,2,3)
%     plot(f1,[x(t),y(t)],z(t)),xlabel('x'),ylabel('y'), title('Gaussian fit')
%     zpred = gaussFun(f1.a,f1.x0,f1.y0,f1.sx,f1.xy,x,y)
    
    % Plot 4: Try with toolbox found online
    [fitresult, zfit, fiterr, zerr, resnorm, rr] = fitGaussian1(x,y,z);
    fitresult %[amp, sx, sxy, sy, xo, yo, zo];;
    zpred = gaussFun(fitresult(1),fitresult(5),fitresult(6),fitresult(2),fitresult(4),x,y);
    subplot(2,2,4)
    plot3(x,y,zpred,'r.'),xlabel('x'),ylabel('y'), title('Fit with Custom func')


end