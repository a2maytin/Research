clc
clear 
close all
Dsim = 0.2;
params.boundary = false;

edges = (0:0.01:.8);
pos = (0.005:0.01:0.795);

r_sim = ssd(Dsim,1,params);
counts_exp = histcounts(r_sim,edges);

y = (0:0.001:.8);
counts_mod = counts_model(Dsim,1,y);



figure
b1 = bar(pos,counts_exp);
hold on
p1 = plot(y,counts_mod,'LineWidth',2);
grid on;

%xlim([0, 0.8]);
xlabel('Single Step Displacement (um)', 'FontSize', 14);
ylabel('Bin Count (Normalized)', 'FontSize', 14);
title('Analytical Formula fit to Simulated Brownian Diffusion','FontSize', 14);
l1 = legend('Simulation','Analytical');
l1.FontSize = 16;

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function counts = counts_model(D, f, y)
dr = 0.01;
tau = 0.021742;
func_array=(1./D')*dr*y/(2*tau).*exp(1./(D'*4*tau)*-y.^2);
counts = f*func_array*44443;
end