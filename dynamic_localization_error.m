clear
% Import the measured trajectories
file = load('SK249-rif_tracksFinal.mat');
traj = file.tracksFinal;
coord = {traj.tracksCoordAmpCG};

pos = {traj.tracksCoordXY};

% convert units to SI units
pixelSize = 160e-9; 
timeStep = 21.742e-3;

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

