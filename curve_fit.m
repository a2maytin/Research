file = load('SK249-rif_tracksFinal.mat');
traj = file.tracksFinal;
coord = {traj.tracksCoordAmpCG};
pos = {traj.tracksCoordXY};
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

limit = .85;
dr = 0.01;
edges = (0:dr:limit);
counts_exp = histcounts(exp, edges);

ydata = counts_exp;

xdata = [dr/2:dr:limit-dr/2];


x0 = [0.069,0.157,0.258,0.472,1.667,     0.263,0.369,0.229,0.111];

x1 = [0.0112966928406229,0.0797728931996451,0.215818583774448,0.558388282946502,2.28165015943585,0.00369426844645392,0.368284086609311,0.503660092603342,0.104579559808946];
counts_test = counts_model(x0, xdata);

x = lsqcurvefit(@counts_model,x0,xdata,ydata);

chi2 = chi_squared(counts_exp,counts_test)

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function counts = counts_model(x,xdata)
D = x(1:5);
f = [x(6:9), 1-sum(x(6:9))];
tau = 0.021742;
dr = 0.01;
y = xdata;
func_array=(1./D')*dr*y/(2*tau).*exp(1./(D'*4*tau)*-y.^2);
counts = f*func_array*44443;
end
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

chisq = sum(pulls.*pulls); %chisq, not per dof
%chisq = chisq/(length(exp)-7); 
end