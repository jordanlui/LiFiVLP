function result = lookupFingerprint(testData,p,k)
    % Fingerprint analysis style lookup
    % Input is query struct, with known points
    % P are the params of a system, 4 Gaussian LED transmitters
    % k are the indices we will query
    % Note error result is absolute (squared)
    
    if nargin <3
        % Then query all points
        k = 1:length(testData.x);
    end
    
    % Simulate param space of Gaussians
    % Create mesh grids
    x = 1:640;
    y = 1:480;
    [x,y] = meshgrid(x,y);
    powerVec = zeros(640*480,4);
    % Calculate power values for each channel
    % This generates the Gaussian database based on parameters
    for i=1:4
        a=p(i,1); x0=p(i,2); y0=p(i,3); sx=p(i,4); sy=p(i,5);
        % This if loop will handle case of threshold or no threshold 
        % (the flat top gaussian threshold)
        if size(p,2)>5
            threshold = p(i,6);
        else
            threshold = 0 ;
        end
        % Calculate power values at each pixel location
        power{i} = gaussFun(a,x0,y0,sx,sy,x,y,threshold);
        % Reshape into a 307200 * 4 matrix for easy power lookups
        powerVec(:,i) = reshape(power{i},[640*480,1]);
    end
    
    % Loop through a different recording and test our accuracy
    M = length(k);
    error = zeros(M,1);
    % Pre-allocate a struct for prediction values
    guess = struct('x',zeros(M,1),'y',zeros(M,1),'d',zeros(M,1));

    % Grab test data. Grab the specific points we will query
    real.x = testData.x(k);
    real.y = testData.y(k);
    real.signal = testData.signal(k,:);
    real.d = sqrt(real.x.^2 + real.y.^2); % Distance from origin. Arbitrary anchor point for our system
    
    %     Run through query points
%     tic
    for i = 1:M

        
        % Method 1 - Localize point based on clustered pick of minimum points 
%         [h,v] = localizeClusteredMinimum(powerVec,real.signal(i,:),10,size(x));
        
        % Method 2 - localize point based on sum square differences
        [h,v] = localizeSumSquareDifference(powerVec,real.signal(i,:),10,size(x));
        
        % Write result into guess vector
        guess.x(i) = h;
        guess.y(i) = v;  

        % Score error as difference in real and guessed position, pixel
        % distance
        error(i) = sqrt((real.x(i) - guess.x(i)).^2 + (real.y(i) - guess.y(i)).^2); 
        
    end
%     toc
    guess.d = sqrt([guess.x].^2 + [guess.y].^2); % Calculate distance values
    
    % Calculate R2 fit based on distance values
    ssRes = sum((real.d - guess.d).^2);
    ssTot = sum((real.d - mean(real.d)).^2);
    R2 = 1-ssRes/ssTot;
    
    result.error = error;
    result.R2 = R2;
    result.guess = guess;
    result.real = real;
    
%     ssTot = ;
end