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
file = load('SK249_tracksFinal.mat');
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

% 4 states
Dbest = [0.030922	0.095971	0.25078	 1.1491]; %vbSPT			
fbest = [0.10211	0.43707	0.39841	0.062407];
Dbest = [0.041,0.131,0.287,0.731];	%limit 0.5, 100 tries -->> 2.1059 +/-0.3807(/100)
fbest = [0.102,0.532,0.272,0.094];
Dbest = [0.041,0.130,0.275,0.680];	%limit 0.8, error min(1,sqrt(exp))
fbest = [0.101,0.523,0.268,0.108];  %------>> 1.8244 +/- 0.2718(/100)

dl1 = 0.001;
dl = 0;
D1 = (Dbest(1)-dl1-dl1-dl1:dl1:Dbest(1)-dl1);
D2 = (Dbest(2)-dl:dl1:Dbest(2)+dl);
D3 = (Dbest(3)-dl:dl1:Dbest(3)+dl);
D4 = (Dbest(4)-dl:dl1:Dbest(4)+dl);
f1 = (fbest(1)-dl:dl1:fbest(1)+dl);
f2 = (fbest(2)-dl:dl1:fbest(2)+dl);
f3 = (fbest(3)-dl:dl1:fbest(3)+dl);
%%
close all

params.boundary = false;


limit = .8;
dr = 0.01;
edges = (0:dr:limit);
counts_exp = histcounts(exp, edges);

% params.boundary = true;
% params.loc = true;
% params.avg = true;


chisq_array = zeros([length(D1),length(D2),length(D3),length(D4),length(f1),length(f2),length(f3)]);
std_array = zeros([length(D1),length(D2),length(D3),length(D4),length(f1),length(f2),length(f3)]);
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


for n = 1:length(D4)
for m = 1:length(D3)
for l = 1:length(D2)
for k = 1:length(D1)
for c = 1:length(f3)
for b = 1:length(f2)
for a = 1:length(f1)
clc
nsim=1; %factor to make simulated data more precise
tries = 100;
f4 = 1-f1(a)-f2(b)-f3(c);
fprintf('\nFitting distribution for D=[%.3f,%.3f,%.3f,%.3f], f=[%.3f,%.3f,%.3f,%.3f]', ...
        D1(k),D2(l),D3(m),D4(n),f1(a),f2(b),f3(c),f4);
fprintf('\nAveraging goodness of fit between exp data and %d simulations, \nCurrently on simulation:  ',tries);

gof2array = zeros([1, tries]);
for j = 1:tries

    Dsim=[D1(k),D2(l),D3(m),D4(n)];
    fsim=[f1(a),f2(b),f3(c),f4];
    sim4 = ssd(Dsim, nsim,params);
    counts_sim4 = 0;
    for i = 1:4
        binlim = nsim*round(17974*fsim(i));
        counti = 1/nsim*histcounts(sim4(i,1:binlim), edges);
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

figure
b1 = bar(counts_exp, 'FaceAlpha', 0.5);
hold on
counts_best = counts_model(Dbest,fbest,limit)*17974;
b4 = bar(counts_best, 'FaceAlpha', 0.5);
title('exp vs model')

figure
b1 = bar(counts_exp, 'FaceAlpha', 0.5);
hold on
b2 = bar(counts_sim4, 'FaceAlpha', 0.5);
title('exp vs sim')

figure
b3 = bar((counts_exp-counts_best)./sqrt(counts_exp));
title('model residuals');

gofmodel = chi_squared(counts_exp,counts_best);

figure
b4 = bar((counts_exp-counts_sim4)./sqrt(counts_exp));
title('simulation residuals');



%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function chisq = chi_squared(exp,sim)  
% calculates chi squared fit to experimental data
residuals = exp-sim;
% assume Poisson statistics
errors = max(1,sqrt(exp));
%errors = 20*ones(1,length(exp));
pulls = residuals./errors;
%figure
%h3 = histogram(pulls)

chisq = sum(pulls.*pulls);
chisq = chisq/(length(exp)-7); %chisq per dof
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