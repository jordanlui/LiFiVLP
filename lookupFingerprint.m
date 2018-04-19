function result = lookupFingerprint(input,p,k)
    % Fingerprint analysis style lookup
    % Input is query struct, with known points
    % P are the params of a system, 4 Gaussian LED transmitters
    % k are the indices we will query
    
    if nargin <3
        % Then query all points
        k = 1:length(input.x);
    end
    
    % Simulate param space of Gaussians
    x = 1:640;
    y = 1:480;
    [x,y] = meshgrid(x,y);
    powerVec = zeros(640*480,4);
    for i=1:4
        a=p(i,1); x0=p(i,2); y0=p(i,3); sx=p(i,4); sy=p(i,5);
        power{i} = gaussFun(a,x0,y0,sx,sy,x,y);
        powerVec(:,i) = reshape(power{i},[640*480,1]);
    end
    
    % Loop through a different recording and test our accuracy
    M = length(k);
    error = zeros(M,1);
    guess = struct('x',zeros(M,1),'y',zeros(M,1));
%     Run through query points
    for i = 1:M
        q.c = [input.signal(k(i),:)];
        q.x = input.x(k(i));
        q.y = input.y(k(i));

        % Look up based on our input values. Select several so we can average
        aSearch = min(powerVec - q.c,[],2);
        [B,ind] = mink(abs(aSearch),10);

        % Convert back to subscript index
        [v,h] = ind2sub(size(x),ind);
        % Average the position
        v = round(mean(v));
        h = round(mean(h));
        guess(i).x = h;
        guess(i).y = v;

        % Evaluate error
        error(i) = sqrt((q.x - guess(i).x).^2 + (q.y - guess(i).y).^2);
    end

    result.error = error;
%     result.guess = guess;
end