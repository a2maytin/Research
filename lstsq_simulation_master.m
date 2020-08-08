% Compute probability distribution of single step displacements from simulated 
% RNaseE trajectories, then perform least squares fit to experimental data
% Find optimal values of fit parameters D and f.

% This code generates a table of chi-squared values for a variety of input 
% parameters but unlike the lstsq_brownian version it does not iterate to 
% find progressively better parameters (doing so would take too much time).
% The code also outputs the standard deviation for each chi-squared value.
% Due to poisson counting error, a simulation with a given set of
% parameters sometimes does better, sometimes worse job of fitting the 
% experimental data.

% Import the measured trajectories
file = load('SK249-rif_tracksFinal.mat'); %SK249-rif.tracksFinal.mat
traj = file.tracksFinal;
coord = {traj.tracksCoordAmpCG};
pos = {traj.tracksCoordXY};

% convert units to SI units,  specific to 249/249-rif data for this project
pixelSize = .16; % in micrometers
timeStep = 21.742e-3;

exp=ones(1,(sum(cellfun(@length, pos))-numel(pos))); %array will store experimental displacements

k = 1; %counter for pooling displacements

for i = 1:numel(pos) %loop through every track
    for j = 1:(length(pos{i})-1) %loop through every displacement
        r_new = sqrt((pos{i}(j+1,1)-pos{i}(j,1))^2+(pos{i}(j+1,2)-pos{i}(j,2))^2);
        exp(k) = r_new*pixelSize;
        k=k+1;
    end
end
%% INITIAL GUESS, SET FIXED HISTOGRAM/SIMULATION PARAMETERS

numStates = 3; %Number of states

limit = .85; % Set the limit for binning to avoid zero count bins
dr = 0.01; % Set the bin width

nsim = 100; % number of simulations to run for each set of parameters

% To aid with initial guessing, the effect of the boundary, localization 
% error and averaging was found to lower the apparent D of brownian
% diffusing particle. D_to_Dsim is designed to take the best brownian fit 
% parameters as input and calculate what would be a good guess for
% simulated parameters.
Ds = [0.088,	0.271, 1.35];
D0 = D_to_Dsim([0.088,	0.271, 1.35]);
f0 = [.471, .486, .043];

% Simulation parameters: turn on boundary, localization error, and 
% averaging due to finite camera exposure
params.boundary = true;
params.loc = true;
params.avg = true;
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
D1 = (.1:.01:.1);
D2 = (.4 :.01 :.4);
D3 = (2.4 :.01 : 2.4);
D4 = 0;
f1 = (.47 :.01 :.47);
f2 = (.49 :.01 :.49);
f3 = 0;


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Dbest = D0;
fbest = f0;
edges = (0:dr:limit);
counts_exp = histcounts(exp, edges);
N = length(exp); 
%%
chisq_array = zeros([length(D1),length(D2),length(D3),length(D4),length(f1),length(f2),length(f3)]);
std_array = zeros([length(D1),length(D2),length(D3),length(D4),length(f1),length(f2),length(f3)]);

for n = 1:length(D4)
for m = 1:length(D3)
for l = 1:length(D2)
for k = 1:length(D1)
for c = 1:length(f3)
for b = 1:length(f2)
for a = 1:length(f1)
clc
Dsim=[D1(k),D2(l),D3(m),D4(n)];
fsim = [f1(a),f2(c),f3(c),0];
fsim(numStates) = 1 - sum(fsim);

fprintf('\nFitting distribution for D=[%.2f,%.2f,%.2f,%.2f], f=[%.2f,%.2f,%.2f,%.2f]', ...
        Dsim(1),Dsim(2),Dsim(3),Dsim(4),fsim(1),fsim(2),fsim(3),fsim(4));
fprintf('\nAveraging goodness of fit between exp data and %d simulations, \nCurrently on simulation:  ',nsim);
gof2array = zeros([1, nsim]);
for j = 1:nsim
    sim4 = ssd(Dsim, N,params);
    counts_sim4 = 0;
    for i = 1:4
        binlim = round(N*fsim(i));
        counti = histcounts(sim4(i,1:binlim), edges);
        counts_sim4 = counts_sim4 + counti;     
    end

    gof2array(j) = chi_squared(counts_exp, counts_sim4);

    if j>10
    fprintf('\b');
    end
    fprintf('\b%d', j);
end
chisq_array(k,l,m,n,a,b,c) = mean(gof2array);
std_array(k,l,m,n,a,b,c) = std(gof2array);
fprintf('  ... Done.\n');
end
end
end
end
end
end
end


disp(chisq_array);
disp(std_array);
chisqPerDof = chisq_array/(length(counts_exp)-(2*numStates-1));
stdPerDof = std_array/(length(counts_exp)-(2*numStates-1));
disp(chisqPerDof);
disp(stdPerDof);

close all

figure
b1 = bar(counts_exp, 'FaceAlpha', 0.5);
hold on
counts_best = counts_model(Ds,f0,timeStep, dr, limit)*N;
b4 = bar(counts_best, 'FaceAlpha', 0.5);
title('Experiment vs Brownian Model')
xlabel('Single Step Displacement (um)')
ylabel('Counts')

figure
b3 = bar((counts_exp-counts_best)./max(1,sqrt(counts_exp)));
title('Brownian Model Residuals');
xlabel('Single Step Displacement (um)')
ylabel('Normalized Residuals')

gofmodel = chi_squared(counts_exp,counts_best);

figure
b1 = bar(counts_exp, 'FaceAlpha', 0.5);
hold on
b2 = bar(counts_sim4, 'FaceAlpha', 0.5);
title('Experiment vs Simulation')
xlabel('Single Step Displacement (um)')
ylabel('Counts')

figure
b4 = bar((counts_exp-counts_sim4)./max(1,sqrt(counts_exp)));
title('Simulation Residuals');
xlabel('Single Step Displacement (um)')
ylabel('Normalized Residuals')

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function chisq = chi_squared(exp,sim)  
residuals = exp-sim;
% assume Poisson statistics
errors = max(1,sqrt(exp));
pulls = residuals./errors;
chisq = sum(pulls.*pulls);
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function counts = counts_model(D, f, tau, dr, limit)
y = [dr/2:dr:limit-dr/2];
func_array=(1./D')*dr*y/(2*tau).*exp(1./(D'*4*tau)*-y.^2);
counts = f*func_array;
end