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
%%
close all

params.boundary = false;
% params.avg = true;
% params.loc = true;
% params.boundary = true;

limit = .8;
dr = 0.01;
edges = (0:dr:limit);
counts_exp = histcounts(exp, edges);

Dbest = [0.088133,	0.24593,	1.3271];
%Dbest = [0.088133,	0.24593,	1.0];
fbest = [0.41606,	0.52508,	0.058861];


%Dbound = D_to_Dsim(Dbest);
%Dbound = [0.02, 0.28, 2.4];
Dbest = [0.08,	0.248,	1.12];
fbound = fbest;
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% D1 = (0.07:0.01:0.08);
% D2 = (0.19:0.01:0.20);
% D3 = (0.42:0.01:0.43);
% D4 = (1.64:0.01:1.65);
% f1 = (0.30:0.1:0.31);
% f2 = (0.53:0.1:0.53);
% f3 = (0.12:0.1:0.12);
% %f4 = (0.04:0.1:0.04);

% D1 = Dbest(1);
% D2 = Dbest(2);
% D3 = Dbest(3);
% D4 = Dbest(4);
% f1 = fbest(1);
% f2 = fbest(2);
% f3 = fbest(3);

D1 = Dbound(1);
D2 = Dbound(2);
D3 = Dbound(3);

f1 = fbound(1);
f2 = fbound(2);


chisq_array = zeros([length(D1),length(D2),length(D3),length(f1),length(f2)]);
std_array = zeros([length(D1),length(D2),length(D3),length(f1),length(f2)]);
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



for m = 1:length(D3)
for l = 1:length(D2)
for k = 1:length(D1)
for b = 1:length(f2)
for a = 1:length(f1)
clc
nsim=1; %factor to make simulated data more precise
tries = 1;
f3 = 1-f1(a)-f2(b);
fprintf('\nFitting distribution for D=[%.2f,%.2f,%.2f], f=[%.2f,%.2f,%.2f]', ...
        D1(k),D2(l),D3(m),f1(a),f2(b),f3);
fprintf('\nAveraging goodness of fit between exp data and %d simulations, \nCurrently on simulation:  ',tries);

gof2array = zeros([1, tries]);
for j = 1:tries

    Dsim=[D1(k),D2(l),D3(m)];
    fsim=[f1(a),f2(b),f3];
    sim3 = ssd(Dsim, nsim,params);
    counts_sim3 = 0;
    for i = 1:3
        binlim = nsim*round(44443*fsim(i));
        counti = 1/nsim*histcounts(sim3(i,1:binlim), edges);
        counts_sim3 = counts_sim3 + counti;     
    end

    gof2array(j) = chi_squared(counts_exp, counts_sim3);

    if j>10
    fprintf('\b');
    end
    fprintf('\b%d', j);

end
chisq_array(k,l,m,a,b) = mean(gof2array);
std_array(k,l,m,a,b) = std(gof2array);
fprintf('  ... Done.\n');
end
end
end
end
end


disp(chisq_array);
disp(std_array);

figure
b1 = bar(counts_exp, 'FaceAlpha', 0.5);
hold on
counts_best = counts_model(Dbest,fbest,limit)*44443;
b4 = bar(counts_best, 'FaceAlpha', 0.5);
title('exp vs model')

figure
b1 = bar(counts_exp, 'FaceAlpha', 0.5);
hold on
b2 = bar(counts_sim3, 'FaceAlpha', 0.5);
title('exp vs sim')

figure
b3 = bar((counts_exp-counts_best)./sqrt(counts_exp));
title('model residuals');

gofmodel = chi_squared(counts_exp,counts_best);

figure
b4 = bar((counts_exp-counts_sim3)./sqrt(counts_exp));
title('simulation residuals');



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