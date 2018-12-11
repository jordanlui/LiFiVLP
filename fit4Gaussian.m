function params = fit4Gaussian(data,symmetric,threshold)
    % Fit a Gaussian to each data channel and return parameters as matrix
    % Newly added feature to choose whether we fit a symmetric or
    % thresholded Gaussian
    % Output: a, x0, y0, sx, sy, threshold
    
    if nargin < 2
        symmetric = 0;
        threshold = 0;
    end
    for j = 1:4
        x = data.x;
        y = data.y;
        z = data.signal(:,j);
        k = 100;

        % Fit Model
        [model,x,y,z] = fitGaussian(x,y,z,k,symmetric,threshold);
        param{j} = coeffvalues(model);
    end
    % Reshape paramaters to matrix. Cell array may prove useful in future
    numFeat = size(param{1},2); % Find number of features, which varies on input param
    params = reshape(cell2mat(param),[numFeat,4])';
end