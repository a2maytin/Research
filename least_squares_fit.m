clc
clear 
close all

% Compute probability distribution of single step displacements from simulated 
% RNaseE trajectories, then perform least squares fit to experimental data
% to find optimal value of the fit parameters (D_fast, D_slow, f_slow)

% Judge goodness of fit using reduced chi-squared statistic

% To estimate the uncertainty in the model parameters, examine 3D grid of 
% ??2 values generated from the unconstrained fits that varied all three 
% parameters Dfast, fslow, and Dslow. Judge fits to be qualitatively poor 
% whenever the value of ??2 exceeds 1.5

% Import the measured trajectories
a = load('SK249-rif_tracksFinal.mat');
b = a.tracksFinal;
pos = {b.tracksCoordXY};
r_exp=ones(1,(sum(cellfun(@length, pos))-numel(pos))); %array will store experimental displacements

%counter for pooling displacements
k = 1;

for i = 1:numel(pos) %loop through every track
    for j = 1:(length(pos{i})-1) %loop through every displacement
        r_new = sqrt((pos{i}(j+1,1)-pos{i}(j,1))^2+(pos{i}(j+1,2)-pos{i}(j,2))^2);
        r_exp(k) = r_new;
        k=k+1;
    end
end

figure('Name','Measured Data');
histogram(r_exp)

% Create simulated trajectories
params1.dt=0.1;  
params1.dt_out=1; 
params1.t_fin=10; 
params1.totR=100;
params1.D=.5;
[rne,tt,params]=random_diffusion_3D_SphCyl(params1);
totR = params.totR;
r_model=ones(1,(size(rne,1)-1)*totR);  %array will store simulated displacements

%counter for pooling displacements
h = 1;

for j = 1:totR %loop through every track
    for i = 1:(size(rne,1)-1) %loop through every displacement
        r_new = sqrt((rne(i+1,3*j-2)-rne(i,3*j-2))^2+(rne(i+1,3*j-1)-rne(i,3*j-1))^2);
        r_model(h) = r_new;
        h=h+1;
    end
end

figure('Name','Simulated Data');
histogram(r_model)




