%Compute probability distribution of single step displacements from simulated 
%RNaseE trajectories, then perform least squares fit to experimental data
%to find optimal value of the fit parameters (D_fast, D_slow, f_slow)

%Judge goodness of fit using reduced chi-squared statistic

%To estimate the uncertainty in the model parameters, examine 3D grid of 
%??2 values generated from the unconstrained fits that varied all three 
%parameters Dfast, fslow, and Dslow. Judge fits to be qualitatively poor 
%whenever the value of ??2 exceeds 1.5

% Import the measured trajectories
a = load('SK249-rif_tracksFinal.mat');
b = a.tracksFinal;
pos = {b.tracksCoordXY};

r_exp=ones(1,(sum(cellfun(@length, pos))-numel(pos))); %experimental displacements

%counter for pooling displacements
k = 1;

for i = 1:numel(pos) %loop through every track
    for j = 1:(length(pos{i})-1) %loop through every displacement
        r_new = sqrt((pos{i}(j+1,1)-pos{i}(j,1))^2+(pos{i}(j+1,2)-pos{i}(j,2))^2);
        r_exp(k) = r_new;
        k=k+1;
    end
end

histogram(r_exp)