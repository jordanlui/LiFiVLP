% Analysis on Jul 26 2017 data
% Static recordings performed at various heights

clc
clear 
close all

load('data_jul26.mat')

%% Calculate average values
for i = 1:length(param)
    asum = sum(cat(3,param{i}{:}),3);
    paramAvg{i} = asum/length(param{i});

end
combined = cat(3,paramAvg{:}); % Average values at each height
clear asum

%% Investigate issues with the 40 Hz transmitter (Channel 3)
k = 100;

% Following will compare the gaussians for each channel at diff heights
clc
close all
for j = [1:4]
    figure(j)
    for i = 1:6
        adata = data{i}(1);
        x = adata.x;
        y = adata.y;
        z = adata.signal(:,j);

        [model,x,y,z] = fitGaussian(x,y,z,k,0,0.3);
        coeffvalues(model)
        % Plot problems
        ax{j}(i) = subplot(3,2,i);
        plot3(x,y,z,'r.')
        plot(model,[x,y],z) % Plot the model surface to see error
        title(sprintf('Gaussian for file %i. a=%.2f,[x0,y0]=%.2f,%.2f. sx=%.2f',i,model.a,model.x0,model.y0,model.sx));
    end
    hlink{j} = linkprop([ax{j}(:)],{'CameraPosition','CameraUpVector'});
end


%% Fit a single and check the values
% clc
% close all

% channel = 3;
% for i = 1:6
%     adata = data{i}(1);
%     x = adata.x;
%     y = adata.y;
%     z = adata.signal(:,channel);
% %     [model,x,y,z] = fitGaussianSymmetric(x,y,z,k);
%     [model,x,y,z] = fitGaussian(x,y,z,k,1);
%     model
% end
% coeffvalues(model)

% close all
% channel = 4;
% i = 5;
% adata = data{i}(1);
% x = adata.x;
% y = adata.y;
% z = adata.signal(:,channel);
% [model,x,y,z] = fitGaussian(x,y,z,k,0,1);
% model
% 
% coeffvalues(model)
% 
% figure()
% plot3(x,y,z,'r.')
% plot(model,[x,y],z) % Plot the model surface to see error

%% Plot out Gaussians and check quality

% Plot out the gaussian trends for each file and visually compare
% close all

% for i = 1:length(param)
%    figure(i)
%    for j = 1:3
%        ax{i}(j) = subplot(1,3,j);
%        plot4Gaussian(param{i}{j})
%    end
%    hlink(i) = linkprop([ax{i}(:)],{'CameraPosition','CameraUpVector'});
%     
% end

