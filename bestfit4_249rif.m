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
pos = (0.005:0.01:0.845);
counts_exp = histcounts(exp, edges);

Dbest = [.077 .208 .50 1.8];
fbest = [.35 .51 .12 .02];


y = (0:0.001:.85);
counts_mod  = counts_model(Dbest,fbest,y);
counts_mod1 = counts_model(Dbest(1),fbest(1),y);
counts_mod2 = counts_model(Dbest(2),fbest(2),y);
counts_mod3 = counts_model(Dbest(3),fbest(3),y);
counts_mod4 = counts_model(Dbest(4),fbest(4),y);


close all
figure
b1 = bar(pos,counts_exp,'FaceAlpha',.6);
b1.BarWidth = 1;
hold on;
p1 = plot(y,counts_mod1,'LineWidth',2);
p2 = plot(y,counts_mod2,'LineWidth',2);
p3 = plot(y,counts_mod3,'LineWidth',2);
p4 = plot(y,counts_mod4,'LineWidth',2);
p  = plot(y,counts_mod,'k','LineWidth',2)
grid on;

%xlim([0, 0.8]);
xlabel('Single Step Displacement (um)', 'FontSize', 14);
ylabel('Bin Count (Normalized)', 'FontSize', 14);
title('Best 4-state Fit: Analytical (Brownian)','FontSize', 14);
l1 = legend([b1,p],{'Experiment (N=44443 steps)','Best 4-state Brownian model'});
l1.FontSize = 16;



%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function counts = counts_model(D, f, y)
dr = 0.01;
tau = 0.021742;
func_array=(1./D')*dr*y/(2*tau).*exp(1./(D'*4*tau)*-y.^2);
counts = f*func_array*44443;
end