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
conversion = .16; %each pixel is 160 nm = .16 um

r_exp=ones(1,(sum(cellfun(@length, pos))-numel(pos))); %array will store experimental displacements

k = 1; %counter for pooling displacements

for i = 1:numel(pos) %loop through every track
    for j = 1:(length(pos{i})-1) %loop through every displacement
        r_new = sqrt((pos{i}(j+1,1)-pos{i}(j,1))^2+(pos{i}(j+1,2)-pos{i}(j,2))^2);
        r_exp(k) = r_new*conversion;
        k=k+1;
    end
end

figure;
h1 = histogram(r_exp)


% Create simulated trajectories
Dsim = 0.07;
steps = 10;
msteps = 50;
tausim = 0.02;

params1.dt=tausim/msteps;  
params1.dt_out=tausim; %time between each frame (exposure time) in s
params1.t_fin=tausim*steps; 
params1.totR=6000;
params1.D=Dsim; %in um^s/s
[rne,tt,params]=random_diffusion_3D_SphCyl(params1);
totR = params.totR;
r_model=ones(1,(size(rne,1)-1)*totR);  %array will store simulated displacements

h = 1; %counter for pooling displacements

for j = 1:totR %loop through every track
    for i = 1:(size(rne,1)-1) %loop through every displacement
        r_new = sqrt((rne(i+1,3*j-2)-rne(i,3*j-2))^2+(rne(i+1,3*j-1)-rne(i,3*j-1))^2);
        r_model(h) = r_new;
        h=h+1;
    end
end

hold on
h2 = histogram(r_model)

h1.Normalization = 'probability';
h1.BinWidth = 0.01;
h2.Normalization = 'probability';
h2.BinWidth = 0.01;
grid on;
xlim([0, 0.8]);
xlabel('Single Step Displacement (um)', 'FontSize', 14);
ylabel('Bin Count', 'FontSize', 14);
title('SSD Histogram: Experimental, Simulated with averaging', 'FontSize', 14);


% Analytical form (free diffusion 3D)
hold on
y = 0:0.01:5;
D = [0.077;0.2;0.43;1.65]; %vbSPT diffusion coeffs
f = [0.31,0.53,0.12,0.04]; %vbSPT occupancies
f2 = 1;
dr = 0.01;
tau=params.dt_out;

%func3=dr*(1./Dsim.^(3/2))*y.^2/(2*sqrt(pi)*tau^(3/2)).*exp((1./Dsim)*-y.^2/(4*tau));

funcv=dr*(1./D)*y/(2*tau).*exp((1./D)*-y.^2/(4*tau)); %experimental vbSPT fit

func2=dr*(1./Dsim)*y/(2*tau).*exp((1./Dsim)*-y.^2/(4*tau));

func=f*funcv;
plot(y,func,'LineWidth',1.5)
plot(y,func2,'LineWidth',1.5)

legend({'r_{exp}','r_{model}','vbPST linear combi.','D_{sim}, f=1'},'FontSize',14)

%% look at abs(displacement in x)

% Create simulated trajectories
Dsim = 0.2;
steps = 10;
msteps = 50;
tausim = 0.02;

params1.dt=tausim/msteps;  
params1.dt_out=tausim; %time between each frame (exposure time) in s
params1.t_fin=tausim*steps; 
params1.totR=6000;
params1.D=Dsim; %in um^s/s
[rne,tt,params]=random_diffusion_3D_SphCyl(params1);
totR = params.totR;
r_model=ones(1,(size(rne,1)-1)*totR);  %array will store simulated displacements

h = 1; %counter for pooling displacements

for j = 1:totR %loop through every track
    for i = 1:(size(rne,1)-1) %loop through every x displacement
        r_new = rne(i+1,3*j-2)-rne(i,3*j-2);
        r_model(h) = r_new;
        h=h+1;
    end
end

figure
h2 = histogram(r_model)
h2.Normalization = 'probability';
h2.BinWidth = 0.01;

%plot a gaussian with stdev sqrt(2Dtau)
hold on
y = -0.4:0.01:0.4;
k = sqrt(2 * Dsim * tausim);
func = 0.01*normpdf(y,0,k);
plot(y,func,'LineWidth',1.5)

