function ax = plotSignal(data)
    % Plot signals as scatter plot
    s = 10; % Plot size
    color = {[0,0,1],[0,1,0],[0,0,0],[1,0,0]};
    for i = 1:4
        scatter3(data.x,data.y,data.signal(:,i),s,color{i})
    end
    
end