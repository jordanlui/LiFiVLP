% Loads data out of specific folder and formats them:
% 1 cell element for each folder
% 1 struct within cell element for each trial (file)

path = '../Data/july26/static/';
addpath(path)

%% Load files from each folder (each height)
% go in each folder in 'static'. Grab height based on regex. Or brute force
% it

files = dir(path);
files = files(3:end);

for i = 1:length(files)
    [token, match] = regexp(files(i).name,'static_h([0-9]+)mm','tokens','match');
    results.height(i) = str2num(token{1}{1});
end

results.directory = {files.name};
%% Processing
% Go through each folder and fit gaussian to each file
% Fit it in a struct with a name for each one

% Loop through each directory
for i = 1:length(results.directory)
    % Grab each of the CSV files
    files = dir(strcat(path,'/',results.directory{i},'/*.csv'));
    N = length(files);
    
    % Load each file into a struct, each folder into a parent cell array
    for j = 1 : N
        aFilename = strcat(files(j).folder,'/',files(j).name);
        [data,allData,normedData,combined,normalizationParameters] = loadTrialDataNormalizeSave(files);
        data{i}(j) = vlcLoader(aFilename);
        param{i}{j} = fit4Gaussian(data{i}(j));
    end

end

save('data_jul26.mat','data','param')