function fig = plot4Gaussian(params)
    % Plot 4 gaussians in the expected test space
    % Note expected space is 640*480; for now we plot in a 640 square 
    %
    
    xmax = 640;
    ymax = 480;
    color = {'b','g','k','r'};
    ptSize = 15;
%     gaussFun = @(a,x0,y0,sx,sy,x,y) a .* exp(-( (x-x0).^2./2./sx^2 + (y-y0).^2./2./sy^2 ));
    k = 10; % Plot every 10 points
    x = 1:k:xmax;
    y = x; % Use same space for x for now
    [x,y] = meshgrid(x,y);
    
    hold on
    for i = 1:4
        a=params(i,1); x0=params(i,2); y0=params(i,3); sx=params(i,4); sy=params(i,5);
        if size(params,2)>5
            threshold = params(i,6);
        else
            threshold = 0;
        end
        z = gaussFun(a,x0,y0,sx,sy,x,y,threshold);
        ax(i) = scatter3(reshape(x,[numel(x),1]),reshape(y,[numel(y),1]),reshape(z,[numel(z),1]),15,color{i});
    end 
    xlabel('x'), ylabel('y'), zlabel('z')
    legend('1 khz', '10 khz', '40 khz', '100 khz')
    title('Plot of Fitted Gaussians')
    
%     % Plot the center points (optional)
%     for i = 1:4
%         a=params(i,1); x0=params(i,2); y0=params(i,3); sx=params(i,4); sy=params(i,5);
%         scatter3(x0,y0,a,50,'m')
%     end 
    
    hold off
    

end