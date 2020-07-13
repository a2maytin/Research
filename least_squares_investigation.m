clear 
clc

y = (0:0.1:10);
x = y + normrnd(0,1,[1,101]);
y = y + normrnd(0,1,[1,101]);

gof = chi_squared(x,y); %should be 1

%% Error estimation
clc
gof12 = 0;
gof1m = 0;
gof1m2 = 0;
tries = 10;

limit = 0.3;
dr = 0.01;
edges = [0:dr:limit];

for i = 1:tries
    sim = ssd(0.15);
    simc = ssd(0.15);
    counts_sim = histcounts(sim, edges);
    counts_simc = histcounts(simc, edges);

    %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    binlim = round(44443*.5);
    counts_simh = histcounts(sim(1:binlim), edges);
    counts_simch = histcounts(simc(1:binlim), edges);
    %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    gof12 = gof12 + chi_squared(counts_sim, counts_simc,1); % should be like 2, and it is
    gof1m = gof1m + chi_squared(counts_model(0.15,1,limit)*44443, counts_sim, 1); % avg 1
    gof1m2 = gof1m2+ chi_squared(counts_model(0.15,1,limit)*44443, counts_simh+counts_simch, 1); % avg 1
end

disp(gof12/tries)
disp(gof1m/tries)
disp(gof1m2/tries)
%% 2 states
%^^^^^^^^^^
fslow = 0.82;
Dslow = 0.14;
Dfast = 0.74;

%% 4 states
clc

gof42 = 0;
gof4m = 0;
tries = 10;

limit = .8;
dr = 0.01;
edges = [0:dr:limit];
D=[.077,.2,.43,1.65];
f=[.31,.53,.12,.04];

for i = 1:tries
    sim4 = ssd(D);
    counts_sim4 = 0;
    for i = 1:4
        binlim = round(44443*f(i));
        counti = histcounts(sim4(i,1:binlim), edges);
        counts_sim4 = counts_sim4 + counti;     
    end
    disp(sum(counts_sim4))

    sim4c = ssd(D);
    counts_sim4c = 0;
    for i = 1:4
        binlim = round(44443*f(i));
        counts_sim4c = counts_sim4c + histcounts(sim4c(i,1:binlim), edges);
    end

    gof42 = gof42 + chi_squared(counts_sim4, counts_sim4c,0); %should be 2
    gof4m = gof4m + chi_squared(counts_model(D,f,limit)*sum(counts_sim4), counts_sim4,0); %should be 1
end

disp(gof42/tries)
disp(gof4m/tries)
%%

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function chisq = chi_squared(x,y,nparam)  
% calculates chi squared fit to experimental data
residuals = x-y;
% assume Poisson statistics
errors = max(1,sqrt(x));
%errors = 20*ones(1,length(exp));
pulls = residuals./errors;
%figure
%h3 = histogram(pulls)

chisq = sum(pulls.*pulls);
chisq = chisq/(length(x)-nparam); %chisq per dof
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