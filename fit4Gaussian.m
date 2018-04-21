function params = fit4Gaussian(data)
    % Fit a Gaussian to each data channel and return parameters
    for j = 1:4
        x = data.x;
        y = data.y;
        z = data.signal(:,j);
        k = 100;

        % Fit Model
        [model,x,y,z] = fitGaussian(x,y,z,k);
        param{j} = coeffvalues(model);
    end
    % Reshape paramaters to matrix. Cell array may prove useful in future
    params = reshape(cell2mat(param),[5,4])';
end