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


limit = .45;
dr = 0.01;
edges = [0:dr:limit];
D=[.07657,.20044,.4329,1.6469];

%f=[.31061,.52968,.12276,.036953]; 
f=[.3182,.52968,.12276,.036953]; % make it "worse"

counts_mod4 = counts_model(D, f, limit)*44443;

counts_exp = histcounts(exp, edges);

vbSPT4 = chi_squared(counts_exp, counts_mod4);

D=[0.13669,	0.73565];
f=[0.82105,	0.17895];

counts_mod2 = counts_model(D, f, limit)*44443;

vbSPT2 = chi_squared(counts_exp, counts_mod2);

close all
%figure
%b1 = bar(counts_exp)
%hold on
%b2 = bar(counts_mod4)
%figure
%b3 = bar((counts_exp-counts_mod4)./sqrt(counts_exp))

% figure
% b1 = bar(counts_exp)
% hold on
% b2 = bar(counts_mod2)

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Dslow = (0.11:0.01:0.11);
Dfast = (0.85:0.01:0.95);
fslow = (0.8:0.1:0.8);

chisq_array = zeros([length(Dslow),length(Dfast),length(fslow)]);
std_array = zeros([length(Dslow),length(Dfast),length(fslow)]);
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

for n = 1:length(Dslow)
    for m = 1:length(Dfast)
        for f = 1:length(fslow)
            clc
            nsim=1; %factor to make simulated data more precise
            tries = 100;
            fprintf('\nFitting distribution for Dslow=%.2f, Dfast=%.2f, fslow=%.1f',Dslow(n),Dfast(m),fslow(f));
            fprintf('\nAveraging goodness of fit between exp data and %d simulations, \nCurrently on simulation:  ',tries);
            
            gof2array = zeros([1, tries]);
            for j = 1:tries
                
                Dsim=[Dslow(n),	Dfast(m)];
                fsim=[fslow(f),	1-fslow(f)];
                sim2 = ssd(Dsim, nsim);
                counts_sim2 = 0;
                for i = 1:2
                    binlim = nsim*round(44443*fsim(i));
                    counti = 1/nsim*histcounts(sim2(i,1:binlim), edges);
                    counts_sim2 = counts_sim2 + counti;     
                end

                gof2array(j) = chi_squared(counts_exp, counts_sim2);

                if j>10
                fprintf('\b');
                end
                fprintf('\b%d', j);

            end
            chisq_array(n,m,f) = mean(gof2array);
            std_array(n,m,f) = std(gof2array);
            fprintf('  ... Done.');
        end
    end
end


disp(chisq_array);
disp(std_array);

% figure
% b1 = bar(counts_exp);
% hold on
% b2 = bar(counts_sim2, 'FaceAlpha', 0.5);
% figure
% b3 = bar((counts_exp-counts_sim2)./sqrt(counts_exp));




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