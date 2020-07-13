% max_likelihood_fit
% perform a log-liklihood fit to the experimental data
% compare to least squares fit; they should match.

% find the likelihood of a simulated data set compared to model:

D=[.077,.2,.43,1.65];
f=[.31,.53,.12,.04];
sim4 = ssd(D);
counts_sim4 = 0;
for i = 1:4
    counts_sim4 = counts_sim4 + f(i)*histcounts(sim4(i,:), edges);
end

sim4c = ssd(D);
counts_sim4c = 0;
for i = 1:4
    counts_sim4c = counts_sim4c + f(i)*histcounts(sim4c(i,:), edges);
end

counts_mod = counts_model(D,f)*44443;

loglik = 
lik = 

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function counts = counts_model(D, f)
tau = 0.021742;
%edges = [0:0.016:0.8];
%y = [0.008:0.016:0.792];
dr = 0.01;
limit = 0.3;
y = [dr/2:dr:limit-dr/2];
func_array=(1./D')*dr*y/(2*tau).*exp(1./(D'*4*tau)*-y.^2);
counts = f*func_array;
end