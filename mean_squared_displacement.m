%% Experimental data
clear all
clc
close all

SPACE_UNITS = 'µm';
TIME_UNITS = 's';

% convert units to SI units
pixelSize = 160e-9; 
dT = 21.742e-3;

a = load('SK249-rif_tracksFinal.mat');
b = a.tracksFinal;
pos = {b.tracksCoordXY};

N_PARTICLES = numel(pos);
tracks = cell(N_PARTICLES, 1);

for i = 1 : N_PARTICLES %N_PARTICLES

    % Time
    time = (0 : size(pos{i},1)-1)' * dT;

    % Position
    X = pos{i}.*pixelSize;

    % Store
    tracks{i} = [time X];

end
clear i X time
close all

%ma.plotTracks;
%ma.labelPlotTracks;

%%
ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
ma = ma.addAll(tracks);
ma = ma.computeMSD;

%% Default MSD (weighted EATA MSD)
close all
%figure
%ma.plotMSD;
figure
ma.plotMeanMSD(gca, true)
[fo, gof] = ma.fitMeanMSD(2);
plot(fo)
ma.labelPlotMSD;
xlim([0, 0.06])
ylim([0, 0.5e-13])
legend off
lerror = 1/2*sqrt(fo.p2+fo.p1*dT/3);
fprintf('Estimation of the dynamic localization error:\n')
fprintf('sigma = %.3g\n', lerror);

%% EATA MSD

figure
ma.plotEATAMSD; 
hold on
[fo2, gof2] = ma.fitEAMSD(2);
plot(fo2)
ma.labelPlotMSD;
xlim([0, 0.06])
ylim([0, 0.5e-13])
legend off
lerror = 1/2*sqrt(fo2.p2+fo2.p1*dT/3);
fprintf('Estimation of the dynamic localization error:\n')
fprintf('sigma = %.3g\n', lerror);

%% EA MSD
ma = ma.computeEAMSD;

figure
ma.plotEAMSD; 
hold on
[fo3, gof3] = ma.fitEAMSD(2);
plot(fo3)
ma.labelPlotMSD;
xlim([0, 0.06])
ylim([0, 0.5e-13])
legend off
lerror = 1/2*sqrt(fo3.p2+fo3.p1*dT/3);
fprintf('Estimation of the dynamic localization error:\n')
fprintf('sigma = %.3g\n', lerror);


%% Fit Individual MSD Cirves, take average D
ma = ma.fitMSD;

good_enough_fit = ma.lfit.r2fit > 0.8;
Dval = ma.lfit.a / 2 / ma.n_dim;
Dmean = mean( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;
Dstd  =  std( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;

fprintf('Estimation of the diffusion coefficient from linear fit of the MSD curves:\n')
fprintf('D = %.3g ± %.3g (mean ± std, N = %d)\n', ...
   Dmean, Dstd, sum(good_enough_fit));

%% Histogram of MSD(n=1)
msd =  {b.MSD};
msd_1 = zeros(1,numel(msd));
for i = 1:numel(msd)
    msd_1(i) = msd{i}(1,1);
end

h = histogram(msd_1)
h.BinWidth = 1e-15;

%% Simulated data
clear all
clc
close all

SPACE_UNITS = 'µm';
TIME_UNITS = 's';

% Create simulated trajectories
Dsim = 0.2;
steps = 100;
msteps = 10;
tausim = 0.02;

params1.dt=tausim/msteps;  
params1.dt_out=tausim; %time between each frame (exposure time) in s
params1.t_fin=tausim*steps; 
params1.totR=1000;
params1.D=Dsim; %in um^s/s
params1.boundary = false;
[rne,params]=rand_diff_3D_SphCyl(params1);
N=params.totR; 
dT=params.dt_out;

tracks2 = cell(N, 1);

for i = 1 : N 

    % Time
    time = (0 : size(rne,1)-1)' * dT;

    % Position
    X = [rne(:,3*i-2), rne(:,3*i-1)];

    % Store
    tracks2{i} = [time X];

end
clear i X time
close all
ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
ma = ma.addAll(tracks2);
%ma.plotTracks;
%ma.labelPlotTracks;
ma = ma.computeMSD;
ma.msd;
%figure
%ma.plotMSD;

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
title('MSD Fit: Simulation w/o Boundary (Brownian)')
% 
% ma = ma.fitMSDandrew;
% 
% good_enough_fit = ma.lfit.r2fit > 0.8;
% Dval=ma.lfit.a;
% Dmean = mean( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;
% Dstd  =  std( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;
% 
% fprintf('Estimation of the diffusion coefficient from linear fit of the MSD curves:\n')
% fprintf('D = %.3g ± %.3g (mean ± std, N = %d)\n', ...
%     Dmean, Dstd, sum(good_enough_fit));


%% Example from author
SPACE_UNITS = 'µm';
TIME_UNITS = 's';

N_PARTICLES = 10;
N_TIME_STEPS = 100;
N_DIM = 2; % 2D

% Typical values taken from studies of proteins diffusing in membranes:
% Diffusion coefficient
D  = 1e-3; % µm^2/s
% Time step between acquisition; fast acquisition!
dT = 0.05; % s,

% Area size, just used to disperse particles in 2D. Has no impact on
% analysis.
SIZE = 2; % µm
k = sqrt(2 * D * dT);
tracks = cell(N_PARTICLES, 1);

for i = 1 : N_PARTICLES

    % Time
    time = (0 : N_TIME_STEPS-1)' * dT;

    % Initial position
    X0 = SIZE .* rand(1, N_DIM);

    % Integrate uncorrelated displacement
    dX = k * randn(N_TIME_STEPS, N_DIM);
    dX(1, :) = X0;
    X = cumsum(dX, 1);

    % Store
    tracks{i} = [time X];

end
clear i X dX time X0
ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
ma = ma.addAll(tracks);
ma.plotTracks;
ma.labelPlotTracks;
ma = ma.computeMSD;
ma.msd;
figure
ma.plotMSD;

figure
ma.plotMeanMSD(gca, true)

[fo, gof] = ma.fitMeanMSD;
plot(fo)
ma.labelPlotMSD;
legend off

ma = ma.fitMSD;

good_enough_fit = ma.lfit.r2fit > 0.8;
Dval=ma.lfit.a;
Dmean = mean( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;
Dstd  =  std( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;

fprintf('Estimation of the diffusion coefficient from linear fit of the MSD curves:\n')
fprintf('D = %.3g ± %.3g (mean ± std, N = %d)\n', ...
    Dmean, Dstd, sum(good_enough_fit));