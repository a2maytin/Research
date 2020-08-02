clc
clear 
close all

% Gradient Descent. NOTE: Doesn't work as well as lstsq_brownian5

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

% 5 states vbSPT
Dbest = [0.075614,	0.14404,	0.29053,	0.44049,	1.7086];		
fbest = [0.24529,   0.38601,	0.24044,	0.093964,	0.034298];

%best so far 1.1358
Dseed = [0.0730,0.1760,0.2360,0.5555,2.4935];
fseed = [0.3140,0.2970,0.2615,0.1100,0.0175];
Dfit = Dseed;
ffit = fseed;

idx = [2,2,2,2,2, 2,2,2,2];
%%
close all

%params.boundary = true;

limit = .8;
dr = 0.01;
edges = (0:dr:limit);
counts_exp = histcounts(exp, edges);
dl = 0.0005;
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

%grad descent
% for i = 1:5
% if idx(i) == 1
%     Dfit(i) = Dfit(i) - dl;
% elseif idx(i) == 2
%     Dfit(i) = Dfit(i);    
% else
%     Dfit(i) = Dfit(i) + dl;        
% end
% end
% %^^^^^^^^^^^^^^^^^^^^^^^^
% for i = 1:4
% if idx(i+5) == 1
%     ffit(i) = ffit(i) - dl;
% elseif idx(i+5) == 2
%     ffit(i) = ffit(i);    
% else
%     ffit(i) = ffit(i) + dl;        
% end
% end

%best so far
f5 = 1 - ffit(1) - ffit(2) - ffit(3) - ffit(4);

Dsim=[Dfit(1),Dfit(2),Dfit(3),Dfit(4),Dfit(5)];
fsim=[ffit(1),ffit(2),ffit(3),ffit(4),f5];
counts_best= counts_model(Dsim,fsim,limit)*44443;

chisq_array = zeros([1,3]);
chisq_array(2) = chi_squared(counts_exp, counts_best);

for i = 1:5
    D = Dfit;
    f = ffit;
    %look down
    D(i) = D(i) - dl;
    counts_down = counts_model(D,f,limit)*44443;
    chisq_array(1) = chi_squared(counts_exp, counts_down);
    %look up
    D(i) = D(i) + 2*dl;
    counts_up = counts_model(D,f,limit)*44443;
    chisq_array(3) = chi_squared(counts_exp, counts_up);
    %find direction of descent
    [~, idx(i)] = min(chisq_array);
end

for i = 1:4
    D = Dfit;
    f = ffit;
    %look down
    f(i) = f(i) - dl;
    counts_down = counts_model(D,f,limit)*44443;
    chisq_array(1) = chi_squared(counts_exp, counts_down);
    %look up
    f(i) = f(i) + 2*dl;
    counts_up = counts_model(D,f,limit)*44443;
    chisq_array(3) = chi_squared(counts_exp, counts_up);
    %find direction of descent
    [~, idx(i+5)] = min(chisq_array);
end

gofmodel = chi_squared(counts_exp,counts_best);

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