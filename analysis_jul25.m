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

combined.max = reshape([allData{5,:,:}],[N,4]);
combined.min = reshape([allData{6,:,:}],[N,4]);

%% Gaussian fit
% Fit each channel on each file, and plot to compare
% close all
% gaussFun = @(a,x0,y0,sx,sy,x,y) a .* exp(-( (x-x0).^2./2./sx^2 + (y-y0).^2./2./sy^2 ));
% % Loop through the files and fit each set of Gaussians and plot
% figure();
% for i = 1:N
%     % Find the 4 Gaussian parameters for a file
%     params{i} = fit4Gaussian(data(i));
%     
%     % Plot the 4 gaussians
%     ax(i) = subplot(2,2,i);
%     plot4Gaussian(params{i});
% end
% 
% % Link the axes so we can visually inspect and compare
% hlink = linkprop([ax(:)],{'CameraPosition','CameraUpVector'});

%% Compare parameters for fit Gaussians

load('params_jul25.mat','params','paramAvg')
% Add up and take average
paramAvg = (params{1} + params{2} + params{3} + params{4})./N;

% save('params_jul25.mat','params','paramAvg')


% Confirm table scale based on max params
% distTx = sqrt( (paramAvg(1,2) - paramAvg(2,2))^2 + (paramAvg(1,3) - paramAvg(2,3))^2)
% scale = distTx / 240
scale = 4/3;

%% Predictions with Gaussian power database
% Loop through the 4 files. Fit Gaussian to a file, and test against
% another file. Report back the error

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

result.R2 = cell2mat(result.R2);
result.errorMed = cellfun(@(x) median(x),result.error);
result.errorMedMM = result.errorMed / scale;
% Following plots each cell to the boxplot
col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});

% figure(),hist(result.error)
figure(), boxplot2(result.error), title('Error distribution on 4 files')