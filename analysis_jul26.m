% Analysis on Jul 26 2017 data
% Static recordings performed at various heights

clc
clear 
close all

addpath('../../../Code/Functions'); % Functions library

% Choose file we review
% investigationName = 'data_threshold_jul26';
investigationName = 'data_jul26_normed_shuffled';
% investigationName = 'data_jul26';
% investigationName = 'data_july25';

% Load file
load(investigationName)
% load('param_threshold_jul26')
label.channel = {'1 khz','10 khz','40 khz','100 khz'};
label.heights = [0, 65, 120, 160, 180, 190];

N_heights = 6;
N_files = 3;
k = 250; % Point frequency to read into model
parameters.pixelScale = 4/3; % The spatial ratio. 1.33 pixel/mm.


%% Predictions with Gaussian power database
data = normedData;
% Perform fitting for each height
tic
for height = 1:N_heights
    for afile = 1:N_files
        [i,j] = permutePair(afile,N_files); % Indices we use to train,test
        p = param{height}{j};
        % Use average parameters instead of a single param
        % This seem to cause code to hang.
%         files = 1:3;
%         files(afile) = [];
%         p = zeros(4,6);
%         for k = files
%             p = p + param{height}{k};
%         end
%         p = p / length(files);
        
        % Grab input data from testing file
        testData = data{height}(afile);
        ind = 1:k:length(testData.x);
        
        % Lookup on predicted values
        tic
        aresult = lookupFingerprint(testData,p,ind);
        toc
        error{height}{afile} = aresult.error;
        R2{height}{afile} = aresult.R2;
        realVal{height}{afile} = aresult.real;
        guess{height}{afile} = aresult.guess;
        
        
    end
    
    sprintf('Finished analysis at height: %i',height)

end
toc
clear result
result.error = error;
result.R2 = R2;
result.realVal = realVal;
result.guess = guess;

save(strcat('result_shuffle_',investigationName),'result')
% 

%% Plot and show results
close all
load(strcat('result_shuffle_',investigationName))
figcount = 1;
for i = 1:N_heights
    result.errorMedians{i} = cellfun(@(x) median(x),result.error{i});
    result.errorMediansMM{i} = cellfun(@(x) median(x)/parameters.pixelScale,result.error{i});
    result.errorCombined{i} = vertcat(result.error{i}{:});
    result.errorCombinedMM{i} = vertcat(result.error{i}{:})/ parameters.pixelScale;
    result.errorMedian(i) = median(result.errorCombined{i});
    result.errorMean(i) = mean(result.errorCombined{i});
    
end
result.errorMeanMM = result.errorMean ./ parameters.pixelScale;
result.errorMedianMM = result.errorMedian ./ parameters.pixelScale;
% Make a boxplot
figure(), boxplot2(result.errorCombinedMM),title('Error at different heights')
xticklabels(label.heights)
xlabel('System Height change (mm)')
ylabel('Error (mm)')
fixfig(gcf,0);
figcount = saveFigIncrement(figcount);

figure(),scatter(label.heights,result.errorMedianMM), title('Height and Median error')
ylabel('Median error (mm)'), xlabel('System height reduction (mm)')
fixfig(gcf,0);
figcount = saveFigIncrement(figcount);

sprintf('Mean error overall is %.2f',mean(result.errorMedianMM))
sprintf('Mean error overall is %.2f',mean(result.errorMeanMM))
