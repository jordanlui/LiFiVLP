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

%% Compare parameters for fit Gaussians


% Add up and take average
paramAvg = (params{1} + params{2} + params{3} + params{4})./N;

save('params_jul25.mat','params','paramAvg')

%% Predictions with Gaussian power database
% Grab a file and a param set
p = params{2};
input = data(1);
k = 1:100:length(input.x);

result = lookupFingerprint(input,p,k);
figure(),hist(result.error)
figure()


% Simulate param space
x = 1:640;
y = 1:480;
[x,y] = meshgrid(x,y);
powerVec = zeros(640*480,4);
for i=1:4
    a=p(i,1); x0=p(i,2); y0=p(i,3); sx=p(i,4); sy=p(i,5);
    power{i} = gaussFun(a,x0,y0,sx,sy,x,y);
    powerVec(:,i) = reshape(power{i},[640*480,1]);
end

%frame up example query
k = 200;
q.c = [input.signal(k,:)];
q.x = input.x(k);
q.y = input.y(k);

% Look up based on our input values. Select several so we can average
aSearch = min(powerVec - q.c,[],2);
[B,ind] = mink(abs(aSearch),10);

% Convert back to subscript index
[v,h] = ind2sub(size(x),ind);
% Average the position
v = round(mean(v));
h = round(mean(h));
guess.x = h;
guess.y = v;

% Evaluate error
error = sqrt((q.x - guess.x).^2 + (q.y - guess.y).^2);
