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

% Normalize data values?



%% Gaussian fitStage Data for 

i = 2;
x = data(i).x;
y = data(i).y;
z = data(i).signal(:,3);
k = 100;
% Grd data and interpolate
% Beware fact that the griddata isn't super robust, and can return NaN
[xi,yi,zi,x,y,z] = prepGaussian(x,y,z,k);
grid.xi = xi; grid.yi = yi; grid.zi = zi;
% [xData, yData, zData] = prepareSurfaceData( x,y,z ); % Built in function to clean data

% save('gaussTest.mat','xi','yi','zi','x','y','z')

% Following method works to fit Gaussian, Apr 26, 2018
x0 = [1 80 180 50 50]; % Guess. [a x0 y0 sx sy]
gaussFun = @(a,x0,y0,sx,sy,x,y) a .* exp(-( (x-x0).^2./2./sx^2 + (y-y0).^2./2./sy^2 ));
compareGaussianFit(x,y,z,x0,gaussFun,grid)

%% Loop analysis
gaussFun = @(a,x0,y0,sx,sy,x,y) a .* exp(-( (x-x0).^2./2./sx^2 + (y-y0).^2./2./sy^2 ));
close all
for i = 1:1
    for j = 1:4
        x = data(i).x;
        y = data(i).y;
        z = data(i).signal(:,j);
        k = 100;
        [xi,yi,zi,x,y,z] = prepGaussian(x,y,z,k);
        grid.xi = xi; grid.yi = yi; grid.zi = zi;
        [amp, ind] = max(z);
        x0 = [1 x(ind) y(ind) 50 50]; % Guess. [a x0 y0 sx sy]
        compareGaussianFit(x,y,z,x0,gaussFun,grid) % 
        suptitle(sprintf('%s, Ch%i',files(i).name,j))
    end
end
