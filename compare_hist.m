clc
clear 
close all

% Compute probability distribution of single step displacements from simulated 
% RNaseE trajectories, then perform least squares fit to experimental data
% to find optimal value of the fit parameters (D_fast, D_slow, f_slow)

% Judge goodness of fit using reduced chi-squared statistic
House
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


% Plot simulated "theo" trajectories 1
Dsim = 0.077;
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
Dsim2 = .2;
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

% Plot simulated "theo" trajectories 3
Dsim = 0.43;
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
r_theo3=ones(1,(size(rne,1)-1)*totR);  %array will store simulated displacements

h = 1; %counter for pooling displacements

for j = 1:totR %loop through every track
    for i = 1:(size(rne,1)-1) %loop through every displacement
        r_new = sqrt((rne(i+1,3*j-2)-rne(i,3*j-2))^2+(rne(i+1,3*j-1)-rne(i,3*j-1))^2);
        r_theo3(h) = r_new;
        h=h+1;
    end
end

% Plot simulated "theo" trajectories 4
Dsim = 1.65;
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
r_theo4=ones(1,(size(rne,1)-1)*totR);  %array will store simulated displacements

h = 1; %counter for pooling displacements

for j = 1:totR %loop through every track
    for i = 1:(size(rne,1)-1) %loop through every displacement
        r_new = sqrt((rne(i+1,3*j-2)-rne(i,3*j-2))^2+(rne(i+1,3*j-1)-rne(i,3*j-1))^2);
        r_theo4(h) = r_new;
        h=h+1;
    end
end


%legend({'Measured SSDs','2 state simulation w/ boundary','2 state simulation no boundary'},'FontSize',14)

%figure
edges = [0:0.01:0.8];
cr_theo = histcounts(r_theo, edges, 'Normalization', 'probability');
cr_theo2 = histcounts(r_theo2, edges, 'Normalization', 'probability');
cr_theo3 = histcounts(r_theo3, edges, 'Normalization', 'probability');
cr_theo4 = histcounts(r_theo4, edges, 'Normalization', 'probability');

cr_full = 0.31*cr_theo + 0.53*cr_theo2 + 0.12*cr_theo3 + 0.04*cr_theo4;
grid on

pos = [0.005:0.01:0.795];
hold on
b2 = bar(pos,cr_full,'BarWidth',1,'FaceAlpha',.6)
%b2.FaceColor = [0.9290 0.6940 0.1250]; %yellow
xlim([0, 0.8]);
xlabel('Single Step Displacement (um)', 'FontSize', 14);
ylabel('Bin Count (Normalized)', 'FontSize', 14);
title('Meausred RNaseE Single Step Displacements','FontSize', 14);
legend({'Measured','Simulated: D=[.077,.2,.43,1.65], f=[.31,.53,.12,.04]'},'FontSize',14)
%%
figure
edges = [0:0.01:0.8];

cr_expm = histcounts(r_expm, edges, 'Normalization', 'probability');

cdiff = (cr_full-cr_expm);

grid on
b1 = bar(pos,cdiff,'FaceAlpha',1)
b1.FaceColor = [0.8500 0.3250 0.0980]; %orange

%xlim([0, 0.8]);
xlabel('Bin Number', 'FontSize', 14);
ylabel('Bin Count Residual (Normalized)', 'FontSize', 14);
title('Simulated Data Deviation from Measured Data','FontSize', 14);
%legend({'2 state w/ boundary','2-state no boundary'},'FontSize',14)

%% 4 way comparison
clc
clear 
close all
% Plot simulated "theo" trajectories 1
Dsim = 0.2;  % in um^2/s
steps = 10;
msteps = 100;
tausim = 0.02; %time between each frame (exposure time) in s
totR = 6000;
sigma = 0.04;

params1.dt=tausim/msteps;  
params1.dt_out=tausim; 
params1.t_fin=tausim*steps; 
params1.totR=totR;
params1.D=Dsim; 
params1.sigma = sigma;
params1.avg = false;
params1.loc= false;
[rne,~,~]=rand_diff_3D_SphCyl(params1);
params1.avg = true;
[rne_a,~,~]=rand_diff_3D_SphCyl(params1); % averaging
params1.loc = true;
[rne_la,~,~]=rand_diff_3D_SphCyl(params1); %loc error and averaging
params1.avg = false;
[rne_l,~,~]=rand_diff_3D_SphCyl(params1); % loc error

r_p=ones(1,(size(rne,1)-1)*totR);  %array will store simulated displacements
r_a=ones(1,(size(rne,1)-1)*totR);  
r_la=ones(1,(size(rne,1)-1)*totR);  
r_l=ones(1,(size(rne,1)-1)*totR);  


h = 1; %counter for pooling displacements

for j = 1:totR %loop through every track
    for i = 1:(size(rne,1)-1) %loop through every displacement
        r_new = sqrt((rne(i+1,3*j-2)-rne(i,3*j-2))^2+(rne(i+1,3*j-1)-rne(i,3*j-1))^2);
        r_p(h) = r_new;
        
        r_new1 = sqrt((rne_a(i+1,3*j-2)-rne_a(i,3*j-2))^2+(rne_a(i+1,3*j-1)-rne_a(i,3*j-1))^2);
        r_a(h) = r_new1;
        
        r_new2 = sqrt((rne_la(i+1,3*j-2)-rne_la(i,3*j-2))^2+(rne_la(i+1,3*j-1)-rne_la(i,3*j-1))^2);
        r_la(h) = r_new2;
        
        r_new3 = sqrt((rne_l(i+1,3*j-2)-rne_l(i,3*j-2))^2+(rne_l(i+1,3*j-1)-rne_l(i,3*j-1))^2);
        r_l(h) = r_new3;
        h=h+1;
    end
end

figure;

h2 = histogram(r_a)
h2.Normalization = 'probability';
h2.BinWidth = 0.01;
hold on
h1 = histogram(r_p)
h1.Normalization = 'probability';
h1.BinWidth = 0.01;

%hold on
%h3 = histogram(r_la)
%h3.Normalization = 'probability';
%h3.BinWidth = 0.01;
hold on
h4 = histogram(r_l)
h4.Normalization = 'probability';
h4.BinWidth = 0.01;
grid on;

xlim([0, 0.8]);
xlabel('Single Step Displacement (um)', 'FontSize', 14);
ylabel('Bin Count (Normalized)', 'FontSize', 14);
title('Localization Error, Averaging Effect on Simulated SSD Histogram','FontSize', 14);

%legend('no loc error, no avg','loc error and avg','FontSize',14)

hold on
dr = 0.01;
y = 0:0.001:5;
func=dr*y/(2*tausim*Dsim + 2*sigma^2).*exp(-y.^2/(4*Dsim*tausim + 4*sigma^2)); %loc error
func2=dr*y/(2*tausim*Dsim - 2.2/3*Dsim*tausim).*exp(-y.^2/(4*Dsim*tausim -4.4/3*Dsim*tausim)); %averaging CHECK THIS
func3=dr*y/(2*tausim*Dsim).*exp(1/(4*Dsim*tausim)*-y.^2); %neither
plot(y,func,'LineWidth',2)
plot(y,func2,'LineWidth',2)
plot(y,func3,'LineWidth',2)
legend('no loc error, no avg','avg','loc error','FontSize',14)
%% 1D displacements

conversion = .16; %each pixel is 160 nm = .16 um

oned=ones(1,2*(sum(cellfun(@length, pos))-numel(pos))); %array will store experimental displacements

k = 1; %counter for pooling displacements

for i = 1:numel(pos) %loop through every track
    for j = 1:(length(pos{i})-1) %loop through every displacement
        one_new = pos{i}(j+1,1)-pos{i}(j,1);
        oned(k) = one_new*conversion;
        k=k+1;
        one_new = pos{i}(j+1,2)-pos{i}(j,2);
        oned(k) = one_new*conversion;
        k=k+1;
    end
end

figure;
h1 = histogram(oned)
h1.Normalization = 'probability';
h1.BinWidth = 0.01;
grid on;
xlim([-.5, 0.5]);
xlabel('\Delta x (um)', 'FontSize', 14);
ylabel('Bin Count (Normalized)', 'FontSize', 14);
title('X Displacements, timestep n=1,2,3','FontSize', 14);

oned=ones(1,2*(sum(cellfun(@length, pos))-numel(pos))); %array will store experimental displacements
%%
oned2=ones(1,2*(sum(cellfun(@length, pos))-2*numel(pos))); %array will store experimental displacements
k = 1; %counter for pooling displacements

for i = 1:numel(pos) %loop through every track
    for j = 1:(length(pos{i})-2) %loop through every displacement
        one_new = pos{i}(j+2,1)-pos{i}(j,1);
        oned2(k) = one_new*conversion;
        k=k+1;
        one_new = pos{i}(j+2,2)-pos{i}(j,2);
        oned2(k) = one_new*conversion;
        k=k+1;
    end
end

oned3=ones(1,2*(sum(cellfun(@length, pos))-3*numel(pos))); %array will store experimental displacements
k = 1; %counter for pooling displacements

for i = 1:numel(pos) %loop through every track
    for j = 1:(length(pos{i})-3) %loop through every displacement
        one_new = pos{i}(j+3,1)-pos{i}(j,1);
        oned3(k) = one_new*conversion;
        k=k+1;
        one_new = pos{i}(j+3,2)-pos{i}(j,2);
        oned3(k) = one_new*conversion;
        k=k+1;
    end
end

hold on
h2 = histogram(oned2)
h2.Normalization = 'probability';
h2.BinWidth = 0.01;
hold on
h3 = histogram(oned3)
h3.Normalization = 'probability';
h2.BinWidth = 0.01;

%% Boundary Comp
clc
clear 
close all
% Plot simulated "theo" trajectories 1
Dsim = 0.2;  % in um^2/s
steps = 10;
msteps = 100;
tausim = 0.02; %time between each frame (exposure time) in s
totR = 6000;
sigma = 0.04;

params1.dt=tausim/msteps;  
params1.dt_out=tausim; 
params1.t_fin=tausim*steps; 
params1.totR=totR;
params1.D=Dsim; 
params1.sigma = sigma;
params1.avg = false;
params1.loc= false;
[rne,~,~]=rand_diff_3D_SphCyl(params1);
params1.boundary = true;
[rne_a,~,~]=rand_diff_3D_SphCyl(params1); % boundary

r_p=ones(1,(size(rne,1)-1)*totR);  %array will store simulated displacements
r_a=ones(1,(size(rne,1)-1)*totR);  


h = 1; %counter for pooling displacements

for j = 1:totR %loop through every track
    for i = 1:(size(rne,1)-1) %loop through every displacement
        r_new = sqrt((rne(i+1,3*j-2)-rne(i,3*j-2))^2+(rne(i+1,3*j-1)-rne(i,3*j-1))^2);
        r_p(h) = r_new;
        
        r_new1 = sqrt((rne_a(i+1,3*j-2)-rne_a(i,3*j-2))^2+(rne_a(i+1,3*j-1)-rne_a(i,3*j-1))^2);
        r_a(h) = r_new1;
       
        h=h+1;
    end
end

figure;


h1 = histogram(r_p)
h1.Normalization = 'probability';
h1.BinWidth = 0.01;
%h1.FaceAlpha = 1;
hold on
h2 = histogram(r_a)
h2.Normalization = 'probability';
h2.BinWidth = 0.01;


grid on;

xlim([0, 0.8]);
xlabel('Single Step Displacement (um)', 'FontSize', 14);
ylabel('Bin Count (Normalized)', 'FontSize', 14);
title('Boundary effect on Simulated SSD Histogram','FontSize', 14);

%legend('no loc error, no avg','loc error and avg','FontSize',14)

%hold on
%dr = 0.01;
%y = 0:0.001:5;
%func3=dr*y/(2*tausim*Dsim).*exp(1/(4*Dsim*tausim)*-y.^2); %neither
%plot(y,func3,'LineWidth',2)
legend('no boundary','boundary','FontSize',14)