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
outName = 'results_threshold_normalized_july25_SSD_Test234.mat';




%% Predictions with Gaussian power database
% Loop through the 4 files. Fit Gaussian to a file, and test against
% another file. Report back the error
% clear
% Original analysis
% load('params_threshold_jul25.mat')
load('data_july25.mat')
load('params_threshold_normalized_jul25') % New parameters with normalized gaussians
path = '../Data/july25/static/';
scale = 4/3;
numTrials = length(data);
testRange = 2:numTrials;
parameters.pixelScale = 4/3; % The spatial ratio. 1.33 pixel/mm.
parameters.numTrials = numTrials;

%% Examine Possible errors in recording 1
% for ind = 1:4
%     
%     figure(2*ind-1)
%     plot(data(ind).time, data(ind).x)
%     figure(2*ind)
%     scatter(data(ind).x, data(ind).y)
% end
%% Visualize Gaussian fit
% close all
% for i = 1:4
% %     Plot the 4 gaussians
%     figure(i);
% %     ax(i) = subplot(2,2,i);
%     plot4Gaussian(params{i});
%     saveas(gcf,sprintf('Figure %i',i));
% end

%% Analysis
% Perform analysis on normalized data instead
dataBackup = data;
data = normedData; 



stepSizeTestData = 250;
for i = 1:length(testRange)
    % Grab a file and a param set, perform a lookup

    % Grab training data from the data struct
%     dataTraining = data(testRange(i));
    % Train on trial1, evaluate on the others
    dataTesting = data(testRange(i));    
    testDataInd = 1:stepSizeTestData:length(dataTesting.x);
    
    % Choose param set to train on
    % Removing test trial, and train on others
    ind = 1:numTrials;
    ind(testRange(i)) = [];
    % Average to combine Gaussian params
    p = zeros(size(params{1}));
    for j = ind
        p = p + params{j};
    end
    p = p ./ size(ind,2);
%     p  = params{1};
    
    % Train/Test: Lookup to find the predicted values. Returns error  R2, guesses
    aresult = lookupFingerprint(dataTesting,p,testDataInd);
    
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
mean(result.errorMedMM)
mean(result.errorMeanMM)
%% Save results

save(outName,'result')

%% Processing on result
clear result
close all
% load('results_july25')
load(outName)

figcount=0;
% Boxplot to compare each trial
col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
figure(), boxplot2(result.errorMM), title('Error distribution'), ylabel('Error (mm)')
fixfig(gcf,0);
figcount = saveFigIncrement(figcount);

% Histogram plot
figure()
for i = 1:length(testRange)
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

figure()
cdfplot(result.errorMM{2})
ylabel('Cumulative distribution'), xlabel('Error (mm)')
fixfig(gcf,0);
figcount = saveFigIncrement(figcount);


%% Experimental
% Look for error trend in distance

figure()
middle = [640/2, 480/2];
for i = 1:length(testRange)
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

% Combine the distance from middle error plot
middle = [640/2, 480/2];

x = cellfun(@(x) x.x, result.real, 'UniformOutput', false)';
x = cat(1,x{:});    
y = cellfun(@(x) x.y, result.real, 'UniformOutput', false)';
y = cat(1,y{:});    
d = ((x - middle(1)).^2 + (y - middle(2)).^2).^(0.5);
e = result.errorMM';
e = cat(1,e{:});
figure()
scatter(d,e)
xlabel('distance from centre (mm)'), ylabel('error (mm)')
suptitle('Error with distance from middle')

fixfig(gcf,0); figcount = saveFigIncrement(figcount);

% Binned heatmap error plot
% figure()
% scatter3(x,y,e);
% Summarize the error data into bins and make a heatmap plot of error
stepSize = 100;
numSteps = round([480/stepSize, 640/stepSize]);
ee = zeros(numSteps);
for i = 1:numSteps(1) % y axis
    for j = 1:numSteps(2) % x axis
        y1 = stepSize*(i-1) + 1;
        y2 = stepSize*i;
        x1 = stepSize*(j-1) + 1;
        x2 = stepSize*j;
        ind = (x1<x & x<x2 & y1<y & y<y2);
        value = nanmean(e(ind));
        if ~isnan(value)
            ee(i,j) = value;
        else
            ee(i,j) = value;
        end
        
    end
end
% Interpolate the missing NaN value
ee(1,6) = mean([ee(1,5),ee(2,5),ee(2,6)]);

figure()
heatmap(gcf,ee,'CellLabelColor','none')
% heatmap(gcf,ee)
title(sprintf('Heatmap of error, mm'))
xlabel('x axis'), ylabel('y axis')
% fixfig(gcf,0); 
figcount = saveFigIncrement(figcount);
clear x y z

clear x y e d ee





% Location of high error results
errorThreshold = 1;
figure()
for i = 1:length(testRange)
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
for i = 1:length(testRange)
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
