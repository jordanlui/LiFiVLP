function ax = plotSignal(data)
    % Plot signals as scatter plot
    s = 10; % Plot size
    color = {[0,0,1],[0,1,0],[0,0,0],[1,0,0]};
    hold on
    
    for i = 1:4
        scatter3(data.x,data.y,data.signal(:,i),s,color{i})
    end
    hold off
    xlabel('x'),ylabel('y'),zlabel('z'),legend('1 khz', '10 khz', '40 khz', '100 khz')
    
end