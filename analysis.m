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
close all
gaussFun = @(a,x0,y0,sx,sy,x,y) a .* exp(-( (x-x0).^2./2./sx^2 + (y-y0).^2./2./sy^2 ));
% Loop through the files and fit each set of Gaussians and plot
figure();
for i = 1:N
    % Find the 4 Gaussian parameters for a file
    params{i} = fit4Gaussian(data(i));
    
    % Plot the 4 gaussians
    ax(i) = subplot(2,2,i);
    plot4Gaussian(params{i});
end

% Link the axes so we can visually inspect and compare
hlink = linkprop([ax(:)],{'CameraPosition','CameraUpVector'});