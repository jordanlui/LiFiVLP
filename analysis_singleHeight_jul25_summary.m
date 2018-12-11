% Goal: Clean up analysis for publication

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
addpath('D:\Users\Jordan\Documents\MATLAB\Library\fmgaussfit')
addpath('D:\Users\Jordan\Documents\MATLAB\Library')
addpath('../../../Code/Functions'); % Functions library
addpath('library')
% addpath('D:\Jordan''s Files\Documents\MATLAB\Libraries\fmgaussfit')
files = dir(strcat(path,'*.csv'));
figcount = 1;
outName = 'results_threshold_paramAvg_july25_SSD.mat';

%% Compare parameters for fit Gaussians

% Add up and take average
% paramAvg = (params{1} + params{2} + params{3} + params{4})./N;
% save('params_threshold_jul25.mat','params','paramAvg')

%% Predictions with Gaussian power database
% Loop through the 4 files. Fit Gaussian to a file, and test against
% another file. Report back the error
% clear
load('params_threshold_jul25.mat')
load('data_july25.mat')
path = '../Data/july25/static/';

scale = 4/3;
numTrials = length(data);

parameters.pixelScale = 4/3; % The spatial ratio. 1.33 pixel/mm.
parameters.numTrials = numTrials;

for i = 1:numTrials
    % Grab a file and a param set, perform a lookup

    % Grab training data from the data struct
    dataTrain = data(i);
    testDataInd = 1:250:length(dataTrain.x);
    
    % Choose param set by removing test trial, and train on others
    ind = 1:numTrials;
    ind(i) = [];
    % Average to combine Gaussian params
    p = zeros(size(params{1}));
    for j = ind
        p = p + params{j};
    end
    p = p ./ size(ind,2);
    
    % Train/Test: Lookup to find the predicted values. Returns error  R2, guesses
    aresult = lookupFingerprint(dataTrain,p,testDataInd);
    
    result.error{i} = aresult.error;
    result.R2{i} = aresult.R2;
    result.real{i} = aresult.real;
    result.guess{i} = aresult.guess;
end

% result.R2 = cell2mat(result.R2);
result.errorMed = cellfun(@(x) median(x),result.error);
result.errorMedMM = result.errorMed / scale;
result.errorMeanMM = cellfun(@(x) mean(x)./scale,result.error);
result.errorMM = cellfun(@(x) x./scale, result.error, 'UniformOutput', false);
result.errorStDev = cellfun(@(x) std(x), result.error, 'UniformOutput', false);

disp('done')
%% Save results

save(outName,'result')

%% Processing on result
clear result
close all
% load('results_july25')
load(outName)


% Boxplot to compare each trial
col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
figure(), boxplot2(result.errorMM), title('Error distribution in 4 trials'), ylabel('Error (mm)')
fixfig(gcf,0);
figcount = saveFigIncrement(figcount);

% Histogram plot
figure()
for i = 1:4
    subplot(2,2,i)
    hist(result.error{i})
end
fixfig(gcf,0);
figcount = saveFigIncrement(figcount);

figure()
hist(result.errorMM{2})
ylabel('Occurrences'), xlabel('Error (mm)')
fixfig(gcf,0);
figcount = saveFigIncrement(figcount);



%% Experimenal
% Look for error trend in distance
close all

figure()
middle = [640/2, 480/2];
for i = 1:4
    subplot(2,2,i)
%     x = result.real{i}.d;
    x = ((result.real{i}.x - middle(1)).^2 + (result.real{i}.y - middle(2)).^2).^0.5;
    y = result.error{i};    
    scatter(x,y)
    xlabel('distance'), ylabel('error')
end
suptitle('Error with distance from middle')
fixfig(gcf,0); figcount = saveFigIncrement(figcount);
clear x y
errorThreshold = 1;

% Location of high error results
figure()
for i = 1:4
    subplot(2,2,i)
    indCheck = result.real{i}.d > errorThreshold; % Grab points above error threshold
    % Grab our plot values
    x = result.real{i}.x(indCheck);
    y = result.real{i}.y(indCheck);
    z = result.error{i}(indCheck);
    scatter(x,y,z) % Plot
    xlim([0,640]), ylim([0,480])
    xlabel('x'), ylabel('y')
end
suptitle(sprintf('Locations of high error points, error > %.2f',errorThreshold))
fixfig(gcf,0); figcount = saveFigIncrement(figcount);
clear x y z

% Single error scatter plot
figure(),i = 1;
    x = result.real{i}.x;
    y = result.real{i}.y;
    z = result.error{i};
    scatter(x,y,z)
    clear x y z
    xlim([0,640]), ylim([0,480])
    xlabel('x'), ylabel('y')
    title('Error as a function of position')
    fixfig(gcf,0); figcount = saveFigIncrement(figcount);
    

% Error with signal power
figure()
for i = 1:4
    subplot(2,2,i)
    indCheck = result.real{i}.d > errorThreshold; % Grab points above error threshold
    % Grab our plot values
%     x = result.real{i}.x(indCheck);
%     y = result.real{i}.y(indCheck);
    e = result.error{i}(indCheck);
    z = sum(result.real{i}.signal,2); % Signal power sum
    z = z(indCheck);
    scatter(z,e);
    ylabel('error (mm)'), xlabel('summed signal power')
end
suptitle(sprintf('Locations of high error points, error > %.2f',errorThreshold))
clear x y z e
fixfig(gcf,0); figcount = saveFigIncrement(figcount);

%% Height analysis jul 26
clear result data
apath = 'result_data_threshold_jul26';
load(apath)
load('data_threshold_jul26')

calcHeightAnalysis(result,label,N_heights)

fixfig(gcf,0);
figcount = saveFigIncrement(figcount);