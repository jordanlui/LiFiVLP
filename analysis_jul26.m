% Analysis on Jul 26 2017 data
% Static recordings performed at various heights

clc
clear 
close all

addpath('../../../Code/Functions'); % Functions library

% Choose file we review
investigationName = 'data_threshold_jul26';
% investigationName = 'data_jul26';
% investigationName = 'data_july25';

% Load file
load(strcat(investigationName,'.mat'))
% load('param_threshold_jul26')
label.channel = {'1 khz','10 khz','40 khz','100 khz'};
label.heights = [0, 65, 120, 160, 180, 190];

N_heights = 6;
N_files = 3;
k = 250; % Point frequency to read into model

%% Predictions with Gaussian power database

% Perform fitting for each height
% for height = 1:N_heights
%     for afile = 1:N_files
%         [i,j] = permutePair(afile,N_files); % Indices we use to train,test
%         p = param{height}{j};
%         input = data{height}(afile);
%         ind = 1:k:length(input.x);
%         
%         % Lookup on predicted values
%         aresult = lookupFingerprint(input,p,ind);
%         result.error{height}{afile} = aresult.error;
%         result.R2{height}{afile} = aresult.R2;
%         result.real{height}{afile} = aresult.real;
%         result.guess{height}{afile} = aresult.guess;
%         
%         
%     end
% 
% end
% 
% save(strcat('result_',investigationName),'result')
% 

%% Plot and show results
close all
load(strcat('result_',investigationName))
for i = 1:N_heights
    result.errorMedians{i} = cellfun(@(x) median(x),result.error{i});
    result.errorCombined{i} = vertcat(result.error{i}{:});
    result.errorMedian{i} = median(result.errorCombined{i});
    result.errorMean{i} = mean(result.errorCombined{i});
end

% Make a boxplot
figure(), boxplot2(result.errorCombined),title('Error at different heights')
xticklabels(label.heights)
xlabel('System Height change (mm)')
ylabel('Error (mm)')
fixfig(gcf,0);

figure(),scatter(label.heights,[result.errorMedian{:}]), title('Height and Median error')
ylabel('Median error (mm)'), xlabel('System height reduction (mm)')
fixfig(gcf,0);

sprintf('Mean error overall is %.2f',mean(cell2mat(result.errorMedian)))
sprintf('Mean error overall is %.2f',mean(cell2mat(result.errorMean)))
