function data = vlcLoader(filename)
    % Load VLC files to a data object
    % Raw signal order is [1,10,40,100khz]
    raw = csvread(filename);
    data.time = raw(:,1);
    data.x = raw(:,2);
    data.y = raw(:,3);
    data.signal = raw(:,4:7);
    data.maxSignal = max(data.signal);
    data.minSignal = min(data.signal);
    data.meanSignal = mean(data.signal);
    data.medianSignal = median(data.signal);
end