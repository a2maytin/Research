% Find the localization error fromn Utrack data structure.
% The value of 24 nm was found which was used in simulation.
% Andrew Maytin

% Import the measured trajectories
file = load('SK249-rif_tracksFinal.mat');
traj = file.tracksFinal;
coord1 = {traj.tracksCoordAmpCG};
pos1 = {traj.tracksCoordXY};
% convert units to SI units
pixelSize = 160e-9; 
timeStep = 21.742e-3;
dT = timeStep;
%%
N_PARTICLES = numel(pos1);
SPACE_UNITS = 'µm';
TIME_UNITS = 's';
tracks = cell(N_PARTICLES, 1);


% Create the tracks to be analyzed by msdanalyzer
for i = 1 : N_PARTICLES %N_PARTICLES
    % Time
    time = (0 : size(pos1{i},1)-1)' * dT;
    % Position
    X = pos1{i}.*pixelSize;
    % Store
    tracks{i} = [time X];
end
clear i X time
close all

% compute the MSD for every track
ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
ma = ma.addAll(tracks);
ma = ma.computeMSD;
ma = ma.fitMSD;
good_enough_fit = ma.lfit.r2fit > 0.8;
Dval = ma.lfit.a / 2 / ma.n_dim;

%% Get Fastest 25% and Slowest 25%
% The cutoffs for the smallest and fastest 25% of particles was found  and
% these were used to find the localization error for the slowest 25% and
% fastest 25% of particles.


% For 249:
% smallCutoff = .13e-13; %smallest 25%
% largeCutoff = 1.7e-13; %largest 25%

% For 249-rif:
smallCutoff = .19e-13; %smallest 25%
largeCutoff = 2e-13; %largest 25%

Dsmallidx = Dval < smallCutoff;
Dbigidx = Dval > largeCutoff;

%%
for i = 1:2
    
if i == 1
    pos = pos1(Dsmallidx);
    coord = coord1(Dsmallidx);
else
    pos = pos1(Dbigidx);
    coord = coord1(Dbigidx);
end
r_expm=ones(1,(sum(cellfun(@length, pos))-numel(pos))); %array will store experimental displacements
stdev_pool=ones(1,2*sum(cellfun(@length, pos))); %array will store experimental displacements

k = 1; %counter for pooling displacements
s = 1; %counter for pooling stdev

for i = 1:numel(pos) %loop through every track
    for j = 1:(length(pos{i})-1) %loop through every displacement
        r_new = sqrt((pos{i}(j+1,1)-pos{i}(j,1))^2+(pos{i}(j+1,2)-pos{i}(j,2))^2);
        r_expm(k) = r_new*pixelSize;
        k=k+1;
    end
end

for i = 1:numel(coord) %loop through every track
    for j = 1:(length(coord{i})/8) %loop through every position
        stdev_new = 0;
        stdev_pool(s) = pixelSize*coord{i}(8*j-3);
        stdev_pool(s+1) = pixelSize*coord{i}(8*j-2);
        s=s+2;
    end
end

mean_stdev = mean(stdev_pool);

disp(mean_stdev*1e9)

% RESULTS:
% 249:
%    23.8616 for slowest 25%
%    24.0239 for fastest 25%

% 249-rif:
%    24.0617 for slowest 25%
%    24.7499 for fastest 25%
end