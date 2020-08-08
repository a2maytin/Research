% The MSD was calculated using the msdanalyzer class which can be found at:
% http://tinevez.github.io/msdanalyzer/
% 
% Since msdanalyzer uses the method of weighted EATAMSD, addditional
% scripts were added to msdanalyzer for non-weighted EATAMSD, as well as
% EAMSD. In the sections below, all three methods are used.
%
% Andrew Maytin
%% Experimental data
% The tracking data is imported
clear 
clc
space_units = 'µm';
time_units = 's';
% convert units to SI units
pixelSize = .16; % in um
dT = 21.742e-3; % in seconds
a = load('SK249-rif_tracksFinal.mat');
b = a.tracksFinal;
pos = {b.tracksCoordXY};
N_PARTICLES = numel(pos);
tracks = cell(N_PARTICLES, 1);
for i = 1 : N_PARTICLES %N_PARTICLES
    time = (0 : size(pos{i},1)-1)' * dT;     % Time
    X = pos{i}.*pixelSize;     % Position
    tracks{i} = [time X];     % Store
end
clear i X time
close all
%% Track length investigation
% The average and standard deviation for the track length was found
tracklength = zeros([1,length(tracks)]);
for i = 1:length(tracks)    
    tracklength(i) = size(tracks{i},1);
end
% for 249:
%track length mean +/- std: 9.47 +/- 8.80
%max track length = 150

% for 249-rif:
%track length mean +/- std: 9.71 +/- 7.71
%max track length = 87
%% Compute MSD Curves
ma = msdanalyzer(2, space_units, time_units);
ma = ma.addAll(tracks);
ma = ma.computeMSD;
%% Default MSD (weighted EATA MSD)
close all
figure
ma.plotMeanMSD(gca, true)
[fo, gof] = ma.fitMeanMSD(2);
plot(fo)
ma.labelPlotMSD;
legend off
D = fo.p1/4;
sigma = 1/2*sqrt(fo.p2+fo.p1*dT/3);
fprintf('Estimation of the dynamic localization error from linear fit of wEATA:\n')
fprintf('sigma = %.3g\n', sigma);
x = sigma^2/(D*dT);
sigma_a = NormMSDInterceptErrorNW(x,2,2122) %not accurate
%%
mode = 4; % Mode 4 is weighted EATAMSD
[fo1, gof1, x] = ma.fitMeanMSDanom(6,fo.p1/4,sigma,mode);
lerror = x(3);
fprintf('Estimation of the dynamic localization error from anomalous fit of wEATA:\n')
fprintf('sigma = %.3g\n', lerror);
disp(x)
t = (0.02:0.02:2);
t1 = 0.0217;
n = round(t/t(1));
r = [x(1), x(2), x(3)];
 y = 2*((2*r(1)*t1^r(2))/((r(2)+2)*(r(2)+1))*((n+1).^(r(2)+2)+(n-1).^(r(2)+2)-2*n.^(r(2)+2)))  ...
     -  (4*r(1)*t1.^r(2)/((r(2)+2)*(r(2)+1)))  +  (2*r(3).^2);
hold on
plot(t,y,'LineWidth',1)
%% non-weighted EATA MSD
figure
ma.plotEATAMSD; 
hold on
[fo2, gof2] = ma.fitEAMSD(2);
plot(fo2)
ma.labelPlotMSD;
legend off
lerror = 1/2*sqrt(fo2.p2+fo2.p1*dT/3);
fprintf('Estimation of the dynamic localization error from linear fit of nwEATA:\n')
fprintf('sigma = %.3g\n', lerror);
%%
mode = 2; % Mode 2 is nonweighted EATAMSD
[fo1, gof1, x] = ma.fitMeanMSDanom(6,fo.p1/4,sigma,mode);
lerror = x(3);
fprintf('Estimation of the dynamic localization error from anomalous fit of wEATA:\n')
fprintf('sigma = %.3g\n', lerror);
disp(x)
t = (0.02:0.02:2);
t1 = 0.0217;
n = round(t/t(1));
r = [x(1), x(2), x(3)];
 y = 2*((2*r(1)*t1^r(2))/((r(2)+2)*(r(2)+1))*((n+1).^(r(2)+2)+(n-1).^(r(2)+2)-2*n.^(r(2)+2)))  ...
     -  (4*r(1)*t1.^r(2)/((r(2)+2)*(r(2)+1)))  +  (2*r(3).^2);
hold on
plot(t,y,'LineWidth',1)
%% EA MSD
ma = ma.computeEAMSD;
figure
ma.plotEAMSD; 
hold on
[fo3, gof3] = ma.fitEAMSD(2);
plot(fo3)
ma.labelPlotMSD;
legend off
lerror = 1/2*sqrt(fo3.p2+fo3.p1*dT/3);
fprintf('Estimation of the dynamic localization error from linear fir of EA:\n')
fprintf('sigma = %.3g\n', lerror);
%%
mode = 3; % Mode 3 is EAMSD
[fo1, gof1, x] = ma.fitMeanMSDanom(6,fo.p1/4,sigma,mode);
lerror = x(3);
fprintf('Estimation of the dynamic localization error from anomalous fit of EA:\n')
fprintf('sigma = %.3g\n', lerror);
disp(x)
t = (0.02:0.02:2);
t1 = 0.0217;
n = round(t/t(1));
r = [x(1), x(2), x(3)];
 y = 2*((2*r(1)*t1^r(2))/((r(2)+2)*(r(2)+1))*((n+1).^(r(2)+2)+(n-1).^(r(2)+2)-2*n.^(r(2)+2)))  ...
     -  (4*r(1)*t1.^r(2)/((r(2)+2)*(r(2)+1)))  +  (2*r(3).^2);
hold on
plot(t,y,'LineWidth',1)
%% Fit Individual MSD Cirves, take average D
ma = ma.fitMSD;
Dval = ma.lfit.a / 2 / ma.n_dim;
%% Fit the Selected D's
clc
clear fo
close all
bounds = [.0,.3]; %Change the bounds to fit different subpopulations
for i = 1:(length(bounds)-1)
% Get subpopulation indices
Didx = (bounds(i)<Dval)&(Dval<bounds(i+1));
% MSD fit the subpopulation using weighted EATAMSD
ma2 = msdanalyzer(2, space_units, time_units);
ma2 = ma2.addAll(tracks(Didx));
ma2 = ma2.computeMSD;
% Linear fit
figure
ma2.plotMeanMSD(gca, true)
[fo, gof] = ma2.fitMeanMSD(2);
plot(fo)
ma2.labelPlotMSD;
legend off
lerror = 1/2*sqrt(fo.p2+fo.p1*dT/3);
bot = bounds(i);
top = bounds(i+1);
fprintf('Fitting Ds %.2f < Dval < %.2f',bot,top)
fprintf('\nsigma = %.3g\n', lerror);
mode = 3;
% Anomalous fit
[fo1, gof1, x] = ma2.fitMeanMSDanom(6,.065,.6,.025,mode);
lerror = x(3);
fprintf('Estimation of the dynamic localization error wEATA:\n')
fprintf('sigma = %.3g\n', lerror);
disp(x)
t = (0.021742:0.021742:2.1742);
t1 = 0.0217;
n = round(t/t(1));
r = x; %Sometimes Lsqnonlin fails to converge
y = 2*((2*r(1)*t1^r(2))/((r(2)+2)*(r(2)+1))*((n+1).^(r(2)+2)+(n-1).^(r(2)+2)-2*n.^(r(2)+2)))  ...
     -  (4*r(1)*t1.^r(2)/((r(2)+2)*(r(2)+1)))  +  (2*r(3).^2);
hold on
plot(t,y,'LineWidth',1)
xlim([0,0.4])
ylim([0,0.2])
end
%% Simulated data

% for SK249: average cell length is 3.7 +/- 0.9 microns, cell width = 1.2 +/- 0.1 microns (mean +/- std). (Note, this is over 91 cells.)
% for SK249+rif: cell Length = 3.0 +/- 0.6 microns, cell width = 1.10 +/- 0.08 microns (mean +/- std). (Note, this is over 42 cells.)
clear all
clc
close all

space_units = 'µm';
time_units = 's';

% Create simulated trajectories
Dsim = 0.2;
steps = 30;
msteps = 100;
tausim = 0.02;
% set the simulation parameters
params1.dt=tausim/msteps;  
params1.dt_out=tausim; %time between each frame (exposure time) in s
params1.t_fin=tausim*steps; 
params1.totR=1000;
params1.D=Dsim; %in um^s/s
params1.boundary = true;
params1.avg = true;
params1.sigma = 0.023;
params1.loc = true;
params1.l0 = 3;
params1.w0 = 1;
[rne,params]=rand_diff_3D_SphCyl(params1);
N=params.totR; 
dT=params.dt_out;

tracks2 = cell(N, 1);

for i = 1 : N 
    time = (0 : size(rne,1)-2)' * dT; %Time
    X = [rne(2:size(rne,1),3*i-2), rne(2:size(rne,1),3*i-1)]; %Position
    tracks2{i} = [time X];   % Store
end
clear i X time
close all
ma = msdanalyzer(2, space_units, time_units);
ma = ma.addAll(tracks2);
ma = ma.computeMSD;
ma.msd;

figure
ma.plotMeanMSD(gca)
[fo, gof] = ma.fitMeanMSD(2);
hold on
plot(fo)

mmsd = ma.getMeanMSD;
t = mmsd(:,1);
x = mmsd(:,2);
dx = mmsd(:,3) ./ sqrt(mmsd(:,4));
errorbar(t, x, dx, 'k')

ma.labelPlotMSD;
legend off
title('MSD Fit: Simulation')
lerror = 1/2*sqrt(fo.p2+fo.p1*dT/3);
fprintf('Estimation of the dynamic localization error:\n')
fprintf('sigma = %.3g\n', lerror);
%% Anomalous fit to simulated data
mode = 4; %weighted EATAMSD
[fo1, gof1, x] = ma.fitMeanMSDanom(6,.2,.8,.023,mode);
lerror = x(3);
fprintf('Estimation of the dynamic localization error wEATA:\n')
fprintf('sigma = %.3g\n', lerror);
disp(x)
t = (0.021742:0.021742:2.1742);
t1 = 0.0217;
n = round(t/t(1));
r = [.11, .8, .023];
 y = 2*((2*r(1)*t1^r(2))/((r(2)+2)*(r(2)+1))*((n+1).^(r(2)+2)+(n-1).^(r(2)+2)-2*n.^(r(2)+2)))  ...
     -  (4*r(1)*t1.^r(2)/((r(2)+2)*(r(2)+1)))  +  (2*r(3).^2);
hold on
plot(t,y,'LineWidth',1)