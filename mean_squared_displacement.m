%% Experimental data
clear all
clc
close all

SPACE_UNITS = 'µm';
TIME_UNITS = 's';

conversion = .16; %each pixel is 160 nm = .16 um

a = load('SK249-rif_tracksFinal.mat');
b = a.tracksFinal;
pos = {b.tracksCoordXY};
dT = 0.021742; % s,

N_PARTICLES = numel(pos);

tracks = cell(N_PARTICLES, 1);

for i = 1 : N_PARTICLES %N_PARTICLES

    % Time
    time = (0 : size(pos{i},1)-1)' * dT;

    % Position
    X = pos{i}.*conversion;

    % Store
    tracks{i} = [time X];

end
clear i X time
close all
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

ma = ma.fitMSDandrew;

good_enough_fit = ma.lfit.r2fit > 0.8;
Dval = ma.lfit.a / 2 / ma.n_dim;
Dmean = mean( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;
Dstd  =  std( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;

fprintf('Estimation of the diffusion coefficient from linear fit of the MSD curves:\n')
fprintf('D = %.3g ± %.3g (mean ± std, N = %d)\n', ...
    Dmean, Dstd, sum(good_enough_fit));

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
[rne,tt,params]=random_diffusion_3D_SphCyl(params1);
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
ma.plotMeanMSDandrew(gca, true)

[fo, gof] = ma.fitMeanMSD;
plot(fo)
ma.labelPlotMSDandrew;
legend off
title('MSD Fit: Simulation w/ localization error and averaging microsteps')

ma = ma.fitMSDandrew;

good_enough_fit = ma.lfit.r2fit > 0.8;
Dval=ma.lfit.a;
Dmean = mean( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;
Dstd  =  std( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;

fprintf('Estimation of the diffusion coefficient from linear fit of the MSD curves:\n')
fprintf('D = %.3g ± %.3g (mean ± std, N = %d)\n', ...
    Dmean, Dstd, sum(good_enough_fit));


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