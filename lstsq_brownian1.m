clc
clear 
close all

% Compute probability distribution of single step displacements from simulated 
% RNaseE trajectories, then perform least squares fit to experimental data
% to find optimal value of the fit parameters D1-4, f1-3

% Judge goodness of fit using reduced chi-squared statistic

% To estimate the uncertainty in the model parameters, examine 7D grid of 
% chi2 values generated from the unconstrained fits that varied all three 
% parameters D1-4, f1-3. Judge fits to be qualitatively poor 
% whenever the value of chi2 exceeds 2

% Import the measured trajectories
file = load('SK249-rif_tracksFinal.mat');
traj = file.tracksFinal;
coord = {traj.tracksCoordAmpCG};

pos = {traj.tracksCoordXY};

% convert units to SI units
%pixelSize = 160e-9; 
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

idx = 1;
%%
close all

%params.boundary = true;

limit = .8;
dr = 0.01;
edges = (0:dr:limit);
counts_exp = histcounts(exp, edges);
dl = 0.005;

% 2 states
Dbest = 0.2;
fbest = 1;

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% BEST: 5.9742
D1b = 0.171;


if idx(1) == 1
    D1b = D1b - dl;
elseif idx(1) == 2
    D1b = D1b;    
else
    D1b = D1b + dl;        
end


%
D1 = (D1b-dl:dl:D1b+dl);


chisq_array = zeros(1,length(D1));
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

fitnum = 0;
z=0;


for k = 1:length(D1)
clc
nsim=1; %factor to make simulated data more precise
tries = 1;
f1 = 1;
fprintf('Fitting distribution for D=[%.2f]', D1(k));

Dsim=[D1(k)];
fsim=[f1];
counts_best = counts_model(Dsim,fsim,limit)*44443;
chisq_array(k) = chi_squared(counts_exp, counts_best);
fitnum = fitnum + 1;
fprintf(repmat('\b',1,z));
msg = num2str(fitnum);
fprintf(msg);
z=numel(msg);
end


fprintf('\n');
% disp(chisq_array);
% 
figure
b1 = bar(counts_exp, 'FaceAlpha', 0.5);
hold on
counts_best = counts_model(Dbest,fbest,limit)*44443;
b4 = bar(counts_best, 'FaceAlpha', 0.5);
title('exp vs model')

figure
b3 = bar((counts_exp-counts_best)./sqrt(counts_exp));
title('model residuals');

gofmodel = chi_squared(counts_exp,counts_best);

[v, linIdx] = min(chisq_array(:));
[idxC{1:ndims(chisq_array)}] = ind2sub(size(chisq_array),linIdx);
idx = linIdx;
disp(v)
disp(idx)
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function chisq = chi_squared(exp,sim)  
% calculates chi squared fit to experimental data
residuals = exp-sim;
% assume Poisson statistics
errors = sqrt(exp);
%errors = 20*ones(1,length(exp));
pulls = residuals./errors;
%figure
%h3 = histogram(pulls)

chisq = sum(pulls.*pulls);
chisq = chisq/(length(exp)-1); %chisq per dof
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function counts = counts_model(D, f, limit)
tau = 0.021742;
%edges = [0:0.016:0.8];
%y = [0.008:0.016:0.792];
dr = 0.01;
y = [dr/2:dr:limit-dr/2];
func_array=(1./D')*dr*y/(2*tau).*exp(1./(D'*4*tau)*-y.^2);
counts = f*func_array;
end