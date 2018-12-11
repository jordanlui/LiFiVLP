% Compare Gaussian distributions from various files
% Worked well on July25 2017 static data.
% Analyzed  Apr 18, 2018
% Future work: Slight tweak to this code could allow easy importing of
% prior data

clc
clear 
close all

path = '../Data/july25/static/';
addpath(path);
addpath('C:\Users\Jordan\Documents\MATLAB\Library\fmgaussfit')
addpath('D:\Jordan''s Files\Documents\MATLAB\Libraries')
% addpath('D:\Jordan''s Files\Documents\MATLAB\Libraries\fmgaussfit')
files = dir(strcat(path,'*.csv'));
%% Load Data
N = 4; % Number of files we want to look at

% Load into structs
for i = 1 : N
    data(i) = vlcLoader(files(i).name);

end

allData = struct2cell(data);
labels = {'1 khz','10 khz','40 khz','100 khz'};
combined.max = reshape([allData{5,:,:}],[N,4]);
combined.min = reshape([allData{6,:,:}],[N,4]);
save('data_july25.mat','data','allData','combined')

%% Gaussian fit on all 4 channels
% Fit each channel on each file, and plot to compare
close all

% Loop through the files and fit each set of Gaussians and plot
figure();
for i = 1:N
    % Find the 4 Gaussian parameters for a file
    params{i} = fit4Gaussian(data(i),0,0.3);
    
    % Plot the 4 gaussians
    ax(i) = subplot(2,2,i);
    plot4Gaussian(params{i});
end

% Link the axes so we can visually inspect and compare
hlink = linkprop([ax(:)],{'CameraPosition','CameraUpVector'});

%% Compare fit with plateau or reg Gaussian
% k = 100;
% clc
% close all
% warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
% for j = 1:4
%     figure(j)
%     labels{j}
%     for i = 1:N
%         adata = data(i);
%         x = adata.x;
%         y = adata.y;
%         z = adata.signal(:,j);
%         [model,x,y,z] = fitGaussian(x,y,z,k,0,0.3);
%         coeffvalues(model)
%         Plot problems
%         ax{j}(i) = subplot(2,2,i);
%         plot3(x,y,z,'r.')
%         plot(model,[x,y],z) % Plot the model surface to see error
%         title(sprintf('Gaussian for %s, file %i. a=%.2f,[x0,y0]=%.2f,%.2f. sx=%.2f',labels{j},i,model.a,model.x0,model.y0,model.sx));
%     end
%     hlink{j} = linkprop([ax{j}(:)],{'CameraPosition','CameraUpVector'});
% end
%% Compare parameters for fit Gaussians


% Add up and take average
paramAvg = (params{1} + params{2} + params{3} + params{4})./N;
save('params_threshold_jul25.mat','params','paramAvg')

% Confirm table scale based on max params
% distTx = sqrt( (paramAvg(1,2) - paramAvg(2,2))^2 + (paramAvg(1,3) - paramAvg(2,3))^2)
% scale = distTx / 240


%% Predictions with Gaussian power database
% Loop through the 4 files. Fit Gaussian to a file, and test against
% another file. Report back the error
% load('params_jul25.mat')
load('params_threshold_jul25.mat')
load('data_july25.mat')
scale = 4/3;
path = '../Data/july25/static/';

for i = 1:4
    % Grab a file and a param set, perform a lookup

    input = data(i);
    k = 1:250:length(input.x);
    if i == 4  % If we reach end of files, loop around
        p = params{1};
    else
        p = params{i+1};
    end

    % Lookup to find the predicted values. Returns error  R2, guesses
    aresult = lookupFingerprint(input,p,k);
    
    result.error{i} = aresult.error;
    result.R2{i} = aresult.R2;
    result.real{i} = aresult.real;
    result.guess{i} = aresult.guess;
end

%% Processing on result
clear
load('results_july25.mat')

result.R2 = cell2mat(result.R2);
result.errorMed = cellfun(@(x) median(x),result.error);
result.errorMedMM = result.errorMed / scale;
% Following plots each cell to the boxplot
col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});

% figure(),hist(result.error)
figure(), boxplot2(result.error), title('Error distribution on 4 files')
% save(strcat(path,'results_threshold_july25.mat'),'result')

% %% Heatmap analysis of error
% % See the spatial distribution of error
% save(strcat(path,'results_july25.mat'),'result')
% figure()
% i = 1;
% plot3(result.real{i}.x,result.real{i}.y,result.error{i},'r.'), title('Spatial Distribution of error')
