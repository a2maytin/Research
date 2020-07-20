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

idx = [2, 2, 2, 2, 2];
%%
close all

%params.boundary = true;

limit = .8;
dr = 0.01;
edges = (0:dr:limit);
counts_exp = histcounts(exp, edges);

dl = 0.001; %change to get better "resolution"

%3 states
Dbest = [0.088133,	0.24593,	1.3271];
fbest = [0.41606,	0.52508,	0.058861];

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%BEST: 
D1b = 0.085;
D2b = 0.264;
D3b = 1.352;
f1b = 0.448;
f2b = 0.508;

if idx(1) == 1
    D1b = D1b - dl;
elseif idx(1) == 2
    D1b = D1b;    
else
    D1b = D1b + dl;        
end
%^^^^^^^^^^^^^^^^^^^^^^^^
if idx(2) == 1
    D2b = D2b - dl;
elseif idx(2) == 2
    D2b = D2b;    
else
    D2b = D2b + dl;        
end
%^^^^^^^^^^^^^^^^^^^^^^^^
if idx(3) == 1
    D3b = D3b - dl;
elseif idx(3) == 2
    D3b = D3b;    
else
    D3b = D3b + dl;        
end
%^^^^^^^^^^^^^^^^^^^^^^^^
if idx(4) == 1
    f1b = f1b - dl;
elseif idx(4) == 2
    f1b = f1b;    
else
    f1b = f1b + dl;        
end
%^^^^^^^^^^^^^^^^^^^^^^^^
if idx(5) == 1
    f2b = f2b - dl;
elseif idx(5) == 2
    f2b = f2b;    
else
    f2b = f2b + dl;        
end

%3^5 = 243

D1 = (D1b-dl:dl:D1b+dl);
D2 = (D2b-dl:dl:D2b+dl);
D3 = (D3b-dl:dl:D3b+dl);
f1 = (f1b-dl:dl:f1b+dl);
f2 = (f2b-dl:dl:f2b+dl);

chisq_array = zeros([length(D1),length(D2),length(D3),length(f1),length(f2)]);
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

fitnum = 0;
z=0;

for m = 1:length(D3)
for l = 1:length(D2)
for k = 1:length(D1)
for b = 1:length(f2)
for a = 1:length(f1)
clc
nsim=1; %factor to make simulated data more precise
tries = 1;
f3 = 1-f1(a)-f2(b);
fprintf('Fitting distribution for D=[%.2f,%.2f,%.2f], f=[%.2f,%.2f,%.2f]\n\n', ...
        D1(k),D2(l),D3(m),f1(a),f2(b),f3);

Dsim=[D1(k),D2(l),D3(m)];
fsim=[f1(a),f2(b),f3];
counts_best = counts_model(Dsim,fsim,limit)*44443;
chisq_array(k,l,m,a,b) = chi_squared(counts_exp, counts_best);
fitnum = fitnum + 1;
fprintf(repmat('\b',1,z));
msg = num2str(fitnum);
fprintf(msg);
z=numel(msg);
end
end
end
end
end


fprintf('\n');
% disp(chisq_array);
% 
% figure
% b1 = bar(counts_exp, 'FaceAlpha', 0.5);
% hold on
% counts_best = counts_model(Dbest,fbest,limit)*44443;
% b4 = bar(counts_best, 'FaceAlpha', 0.5);
% title('exp vs model')
% 
% figure
% b3 = bar((counts_exp-counts_best)./sqrt(counts_exp));
% title('model residuals');

gofmodel = chi_squared(counts_exp,counts_best);

[v, linIdx] = min(chisq_array(:));
[idxC{1:ndims(chisq_array)}] = ind2sub(size(chisq_array),linIdx);
idx = cell2mat(idxC);
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