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

r_expm=ones(1,(sum(cellfun(@length, pos))-numel(pos))); %array will store experimental displacements

k = 1; %counter for pooling displacements

for i = 1:numel(pos) %loop through every track
    for j = 1:(length(pos{i})-1) %loop through every displacement
        r_new = sqrt((pos{i}(j+1,1)-pos{i}(j,1))^2+(pos{i}(j+1,2)-pos{i}(j,2))^2);
        r_expm(k) = r_new*conversion;
        k=k+1;
    end
end

figure;
h1 = histogram(r_expm)
h1.Normalization = 'probability';
h1.BinWidth = 0.01;
grid on;
xlim([0, 0.8]);
xlabel('Single Step Displacement (um)', 'FontSize', 14);
ylabel('Bin Count', 'FontSize', 14);
title('Simulation D=[.214, .75], f=[.8, .2], Boundary vs. No Boundary vs. Measured data','FontSize', 14);

% Plot simulated "exp" trajectories 1
Dsim = 0.14;
steps = 10;
msteps = 50;
tausim = 0.02;

params1.dt=tausim/msteps;  
params1.dt_out=tausim; %time between each frame (exposure time) in s
params1.t_fin=tausim*steps; 
params1.totR=6000;
params1.D=Dsim; %in um^s/s
[rne,tt,params]=rand_diff_3D_SphCyl_exp(params1);
totR = params.totR;
r_exp=ones(1,(size(rne,1)-1)*totR);  %array will store simulated displacements

h = 1; %counter for pooling displacements

for j = 1:totR %loop through every track
    for i = 1:(size(rne,1)-1) %loop through every displacement
        r_new = sqrt((rne(i+1,3*j-2)-rne(i,3*j-2))^2+(rne(i+1,3*j-1)-rne(i,3*j-1))^2);
        r_exp(h) = r_new;
        h=h+1;
    end
end

% Plot simulated "exp" trajectories 2
Dsim2 = 0.75;
steps = 10;
msteps = 50;
tausim = 0.02;

params1.dt=tausim/msteps;  
params1.dt_out=tausim; %time between each frame (exposure time) in s
params1.t_fin=tausim*steps; 
params1.totR=6000;
params1.D=Dsim2; %in um^s/s
[rne,tt,params]=rand_diff_3D_SphCyl_exp(params1);
totR = params.totR;
r_exp2=ones(1,(size(rne,1)-1)*totR);  %array will store simulated displacements

h = 1; %counter for pooling displacements

for j = 1:totR %loop through every track
    for i = 1:(size(rne,1)-1) %loop through every displacement
        r_new = sqrt((rne(i+1,3*j-2)-rne(i,3*j-2))^2+(rne(i+1,3*j-1)-rne(i,3*j-1))^2);
        r_exp2(h) = r_new;
        h=h+1;
    end
end

% Plot simulated "theo" trajectories 1
Dsim = 0.14;
steps = 10;
msteps = 50;
tausim = 0.02;

params1.dt=tausim/msteps;  
params1.dt_out=tausim; %time between each frame (exposure time) in s
params1.t_fin=tausim*steps; 
params1.totR=6000;
params1.D=Dsim; %in um^s/s
[rne,tt,params]=rand_diff_3D_SphCyl_theo(params1);
totR = params.totR;
r_theo=ones(1,(size(rne,1)-1)*totR);  %array will store simulated displacements

h = 1; %counter for pooling displacements

for j = 1:totR %loop through every track
    for i = 1:(size(rne,1)-1) %loop through every displacement
        r_new = sqrt((rne(i+1,3*j-2)-rne(i,3*j-2))^2+(rne(i+1,3*j-1)-rne(i,3*j-1))^2);
        r_theo(h) = r_new;
        h=h+1;
    end
end

% Plot simulated "theo" trajectories 2
%Dsim2 = 0.73565;
Dsim2 = .75;
steps = 10;
msteps = 50;
tausim = 0.02;

params1.dt=tausim/msteps;  
params1.dt_out=tausim; %time between each frame (exposure time) in s
params1.t_fin=tausim*steps; 
params1.totR=6000;
params1.D=Dsim2; %in um^s/s
[rne,tt,params]=rand_diff_3D_SphCyl_theo(params1);
totR = params.totR;
r_theo2=ones(1,(size(rne,1)-1)*totR);  %array will store simulated displacements

h = 1; %counter for pooling displacements

for j = 1:totR %loop through every track
    for i = 1:(size(rne,1)-1) %loop through every displacement
        r_new = sqrt((rne(i+1,3*j-2)-rne(i,3*j-2))^2+(rne(i+1,3*j-1)-rne(i,3*j-1))^2);
        r_theo2(h) = r_new;
        h=h+1;
    end
end

r3 = [r_theo,r_theo,r_theo,r_theo,r_theo,r_theo,r_theo,r_theo,...
      r_theo2,r_theo2]; %ratio: 8 to 2
hold on
h3 = histogram(r3, 'FaceColor', [0.8 0 0])
h3.Normalization = 'probability';
h3.BinWidth = 0.01;

r4 = [r_exp,r_exp,r_exp,r_exp,r_exp,r_exp,r_exp,r_exp,...
      r_exp2,r_exp2]; %ratio: 8 to 2
hold on
h4 = histogram(r4)
h4.Normalization = 'probability';
h4.BinWidth = 0.01;

% Analytical form (free diffusion 3D)
hold on
y = 0:0.001:5;
D = [0.077;0.2;0.43;1.65]; %vbSPT diffusion coeffs
f = [0.31,0.53,0.12,0.04]; %vbSPT occupancies
D2 = [0.13669;.73565]; %vbSPT diffusion coeffs
%D2 = [0.14;.6]; %vbSPT diffusion coeffs
f2 = [0.82105,.17895]; %vbSPT occupancies
%f2 = [0.8,.2]; %test occupancies
dr = 0.01;
tau=params.dt_out;

%func3=dr*(1./Dsim.^(3/2))*y.^2/(2*sqrt(pi)*tau^(3/2)).*exp((1./Dsim)*-y.^2/(4*tau));

funcv=dr*(1./D)*y/(2*tau).*exp((1./D)*-y.^2/(4*tau)); %experimental vbSPT fit
funcv2=dr*(1./D2)*y/(2*tau).*exp((1./D2)*-y.^2/(4*tau)); %worse experimental vbSPT fit

%func2=dr*(1./Dsim)*y/(2*tau).*exp((1./Dsim)*-y.^2/(4*tau));

%hold on
%func=f*funcv;
%plot(y,func,'LineWidth',2)

%hold on
%func2=f2*funcv2;
%plot(y,func2,'LineWidth',2)

legend({'Measured SSDs','2 state simulation w/ boundary','2 state simulation no boundary'},'FontSize',14)

figure
edges = [0:0.01:0.8];
cr_expmnn = histcounts(r_expm, edges);
cr_expm = histcounts(r_expm, edges, 'Normalization', 'probability');
nfactor = cr_expmnn(1)/cr_expm(1);
cr3 = histcounts(r3, edges, 'Normalization', 'probability');
cr4 = histcounts(r4, edges, 'Normalization', 'probability');
cdiff3 = nfactor*(cr3-cr_expm);
cdiff4 = nfactor*(cr4-cr_expm);
grid on
b1 = bar(cdiff3,'FaceAlpha',1)
b1.FaceColor = [0.8500 0.3250 0.0980]; %orange
hold on
b2 = bar(cdiff4,'FaceAlpha',.7)
b2.FaceColor = [0.9290 0.6940 0.1250]; %yellow
%xlim([0, 0.8]);
xlabel('Bin Number', 'FontSize', 14);
ylabel('Residual Count', 'FontSize', 14);
title('Histogram Count Residuals','FontSize', 14);
legend({'2 state w/ boundary','2-state no boundary'},'FontSize',14)
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
r=ones(1,(size(rne,1)-1)*totR);  %array will store simulated displacements

h = 1; %counter for pooling displacements

for j = 1:totR %loop through every track
    for i = 1:(size(rne,1)-1) %loop through every x displacement
        r_new = rne(i+1,3*j-2)-rne(i,3*j-2);
        r(h) = r_new;
        h=h+1;
    end
end

figure
h2 = histogram(r)
h2.Normalization = 'probability';
h2.BinWidth = 0.01;

%plot a gaussian with stdev sqrt(2Dtau)
hold on
y = -0.4:0.01:0.4;
k = sqrt(2 * Dsim * tausim);
func = 0.01*normpdf(y,0,k);
plot(y,func,'LineWidth',1.5)

