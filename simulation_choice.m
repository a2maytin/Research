close all
clc
limit = 0.8;
dr = 0.01;
edges = (0:dr:limit);
%counts_exp = histcounts(exp, edges);


Dm = [.075,.075];
D = D_to_Dsim([.075,.075]);


binlim =0.8;

params.boundary = true;
%params.boundary = false;
params.avg = true;
params.loc = true;

sim_b = ssd(D, 44443, params);

counts_b1 = histcounts(sim_b(1,:), edges);  
counts_b2 = histcounts(sim_b(2,:), edges);  



counts_free = histcounts(ssd(D,nsim),edges);



figure
b1 = bar(counts_b1, 'FaceAlpha', 0.5);
hold on
b2 = bar(counts_b2, 'FaceAlpha', 0.5);
%b3 = bar(counts_free, 'FaceAlpha', 0.5);
%p = plot(counts_model(D,1,limit, 44443*nsim));
counts_mod1 = counts_model(Dm(1),1,limit, 44443*nsim);
counts_mod2 = counts_model(Dm(2),1,limit, 44443*nsim);
p1 = plot(counts_mod1);
p2 = plot(counts_mod2);

figure
b4 = bar((counts_mod1-counts_b1)./sqrt(counts_mod1));
hold on
b5 = bar((counts_mod2-counts_b2)./sqrt(counts_mod2));

disp(chi_squared(counts_mod1,counts_b1))
disp(chi_squared(counts_mod2,counts_b2))

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
function counts = counts_model(D, f, limit, counts)
tau = 0.021742;
%edges = [0:0.016:0.8];
%y = [0.008:0.016:0.792];
dr = 0.01;
scale=counts;
y = [dr/2:dr:limit-dr/2];

func_array=scale*(1./D')*dr*y/(2*tau).*exp(1./(D'*4*tau)*-y.^2);
counts = f*func_array;
end