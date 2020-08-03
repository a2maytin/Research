% Find the localization error fromn Utrack data structure.

% Import the measured trajectories
file = load('SK249_tracksFinal.mat');
traj = file.tracksFinal;
coord = {traj.tracksCoordAmpCG};
pos = {traj.tracksCoordXY};
% convert units to SI units
pixelSize = 160e-9; 
timeStep = 21.742e-3;
dT = timeStep;
%%
N_PARTICLES = numel(pos);
SPACE_UNITS = 'µm';
TIME_UNITS = 's';
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

ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
ma = ma.addAll(tracks);
ma = ma.computeMSD;
ma = ma.fitMSD;
good_enough_fit = ma.lfit.r2fit > 0.8;
Dval = ma.lfit.a / 2 / ma.n_dim;
Dsmallidx = Dval < 0.1e-12;
Dbigidx = Dval > 1e-12;

%%
%^^^^^^^^^^^^^^^
pos = pos(Dsmallidx);
coord = coord(Dsmallidx);
% pos = pos(Dbigidx);
% coord = coord(Dbigidx);
%^^^^^^^^^^^^^^^
%%
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

% mean for D < 0.1e-12 = 2.3442 e-8   vs 5.06 e-8 (4.32 for goodenough r^2 only)
% mean for D > 0.1e-12 = 2.4225 e-8
% mean for all D = 2.3750 e-8        vs   5.07 e-8

%%
Dsmall = 0.1e-12;
sigma_small = 4.32e-8;
x = sigma_small^2/(Dsmall*dT);
sigma_a_small = NormMSDInterceptErrorNW(x,2,100) %what is N (100?)
