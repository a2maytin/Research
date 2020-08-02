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

idx = 2;

close all

%params.boundary = true;

limit = .8;
dr = 0.01;
edges = (0:dr:limit);
counts_exp = histcounts(exp, edges);
dl = 0.005;

% 5 states

Dvbspt = [0.027022	0.076286	0.1796	0.54041	1.8629];
fvbspt = [0.075221	0.29481	0.47998	0.13112	0.018861];
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%best so far 0.9189/65.2419 %good
% D=[0.031,0.091,0.186,0.537,2.484], f=[0.052,0.283,0.484,0.175,0.006]


%limit .8
%D=[0.040,0.103,0.182,0.535,4???], f=[0.093,0.248,0.472,0.182,0.005]
%0.8953/?

%limit .5

%D=[0.037,0.097,0.192,0.561,1.672], f=[0.074,0.279,0.481,0.158,0.008]
%1.0094


%bigg bins 0.015 ---> 1.1540

%limit .5;
%D=[0.037,0.091,0.184,0.526,2.414], f=[0.066,0.250,0.499,0.169,0.016]
%.9389

%limit 1;
%D=[0.037,0.093,0.179,0.531,2.498], f=[0.069,0.245,0.494,0.187,0.005]
%.8 something

%limit .8
%D=[0.037,0.093,0.179,0.536,6.748], f=[0.069,0.245,0.494,0.188,0.004]
%1.0364

Dbest=[0.037,0.093,0.179,0.536,6.748];
fbest=[0.069,0.245,0.494,0.188,0.004];



D1b = Dbest(1);
D2b = Dbest(2);
D3b = Dbest(3);
D4b = Dbest(4);
D5b = Dbest(5);
f1b = fbest(1);
f2b = fbest(2);
f3b = fbest(3);
f4b = fbest(4);

%%
for n = 1:100

%^^^^^^^^^^^^^^^^^^^^^^^^
if idx(1) == 1
    D5b = D5b - 1*dl;
elseif idx(1) == 2
    D5b = D5b;    
else
    D5b = D5b + 1*dl;        
end
%^^^^^^^^^^^^^^^^^^^^^^^^

%3^9 = 19683
D5 = (D5b-dl:dl:D5b+dl);

chisq_array = zeros([1,length(D5)]);
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

fitnum = 0;
z=0;

for o = 1:length(D5)
clc
nsim=1; %factor to make simulated data more precise
tries = 1;
f5 = 1-f1b-f2b-f3b-f4b;
fprintf('Fitting distribution for D=[%.3f,%.3f,%.3f,%.3f,%.3f], f=[%.3f,%.3f,%.3f,%.3f,%.3f]\n\n', ...
        D1b,D2b,D3b,D4b,D5(o),f1b,f2b,f3b,f4b,f5);

Dsim=[D1b,D2b,D3b,D4b,D5(o)];
fsim=[f1b,f2b,f3b,f4b,f5];
counts_best = counts_model(Dsim,fsim,limit)*17974;
chisq_array(o) = chi_squared(counts_exp, counts_best);
fitnum = fitnum + 1;
fprintf(repmat('\b',1,z));
msg = num2str(fitnum);
fprintf(msg);
z=numel(msg);
end

fprintf('\n');
% disp(chisq_array);

[v, linIdx] = min(chisq_array(:));
idx = linIdx;
disp(v)
disp(idx)
f5b = 1- f1b-f2b-f3b-f4b;
fprintf('Best was D=[%.3f,%.3f,%.3f,%.3f,%.3f], f=[%.3f,%.3f,%.3f,%.3f,%.3f]\n\n', ...
        D1b,D2b,D3b,D4b,D5(idx(1)),...
        f1b,f2b,f3b,f4b,f5b);
    
if f5b <= 0
  fprintf('f5 zeroed out.');
  break  
end
    
end

%%
f5 = 1-f1b-f2b-f3b-f4b;
counts_best = counts_model([D1b,D2b,D3b,D4b,D5b],[f1b,f2b,f3b,f4b,f5],limit)*17974;
figure
b1 = bar(counts_exp, 'FaceAlpha', 0.5);
hold on
b4 = bar(counts_best, 'FaceAlpha', 0.5);
title('exp vs model')

figure
b3 = bar((counts_exp-counts_best)./sqrt(counts_exp));
title('model residuals');
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
chisq = chisq/(length(exp)-9); %chisq per dof
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