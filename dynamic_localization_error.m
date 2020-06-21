% To estimate the dynamic localization errors sigma_fast and sigma_slow for
% RNaseE, use six-step trajectories to form the distribution of the mean of 
% six one-step estimates of D

% Import the measured trajectories
a = load('SK249-rif_tracksFinal.mat');
b = a.tracksFinal;
pos = {b.tracksCoordXY};


%% Find the fast and slow cutoff

tau=21.742; %time between each frame (exposure time) in ms

D6_mean=1/(24*tau)*sum()


%% Form estimate of D for slowest 10% of trajectories and the fastest 10% of trajectories


%% Estimate D for single state (249-rif)