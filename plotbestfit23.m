% This script was used to make figures for 


% Import SK249-rif
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
counts_exprif = histcounts(exp, edges)/length(exp);
%% Import SK249
file = load('SK249_tracksFinal.mat');
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
counts_exp249 = histcounts(exp, edges)/length(exp);
%%
% 249-rif
Dbestrif = [.115 .440 1 1 ];
fbestrif = [.707 .293 0 0];
%
%Dbestrif = [.088 .271 1.35 1 ];
%fbestrif = [.471 .486 .043 0];


y = (0:0.001:.85);
%249-rif
counts_rif  = counts_model(Dbestrif,fbestrif,y);
counts_rif1 = counts_model(Dbestrif(1),fbestrif(1),y);
counts_rif2 = counts_model(Dbestrif(2),fbestrif(2),y);
counts_rif3 = counts_model(Dbestrif(3),fbestrif(3),y);
counts_rif4 = counts_model(Dbestrif(4),fbestrif(4),y);



%Top plot
close all
fig = figure;
s1 = subplot(2,1,1);
b1 = bar(pos,counts_exprif,'FaceAlpha',.0);
b1.BarWidth = 1;
hold on;
p1 = plot(y,counts_rif1,'b','LineWidth',2);
hold on
p2 = plot(y,counts_rif2,'b','LineWidth',2);
%p3 = plot(y,counts_rif3,'b','LineWidth',2);
%p4 = plot(y,counts_rif4,'b','LineWidth',2);
prif = plot(y,counts_rif,'k','LineWidth',2)
grid on;
set(gca,'xticklabel',{[]})

xlim([0, 0.5]);
% title('Best Fit: Brownian Model','FontSize', 14);
%l1 = legend([b1,p],{'Experiment (N=44443 steps)','Best 4-state Brownian model'});
% l1 = legend([prif,p1],{'SK249-rif 2-state','subpopulations'});
% l1.FontSize = 16;
hold on
subplot(2,1,2);
% 249
Dbest = [.055 .173 .60 1];
fbest = [.20 .639 .161 0];
%249
counts_249  = counts_model(Dbest,fbest,y);
counts_2491 = counts_model(Dbest(1),fbest(1),y);
counts_2492 = counts_model(Dbest(2),fbest(2),y);
counts_2493 = counts_model(Dbest(3),fbest(3),y);
counts_2494 = counts_model(Dbest(4),fbest(4),y);
b2 = bar(pos,counts_exp249,'FaceAlpha',.0);
b2.BarWidth = 1;
hold on
p5 = plot(y,counts_2491,'r','LineWidth',2);
p6 = plot(y,counts_2492,'r','LineWidth',2);
p7 = plot(y,counts_2493,'r','LineWidth',2);
p8 = plot(y,counts_2494,'r','LineWidth',2);
p  = plot(y,counts_249,'k','LineWidth',2);
% legend
% l1 = legend([p,p5],{'SK249 3-state','subpopulations'});
% l1.FontSize = 16;
grid on;
xlim([0, 0.5]);
a = get(subplot(2,1,2),'position'); % get the current position property
a(2) = 2*a(2) ;      
set(subplot(2,1,2),'position',a);   % set the new position
% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'P(r, t = 21.742 ms)', 'FontSize', 14);
xlabel(han,'Single Step Displacement (um)', 'FontSize', 14);
%title(han,'yourTitle');
y = get(gca,'ylabel'); % handle to the label object
p = get(y,'position'); % get the current position property
p(1) = 1.5*p(1) ;      
set(y,'position',p);   % set the new position
y = get(gca,'xlabel'); % handle to the label object
p = get(y,'position'); % get the current position property
p(2) = -1.5*p(2) ;      
set(y,'position',p);   % set the new position



%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function counts = counts_model(D, f, y)
dr = 0.01;
tau = 0.021742;
func_array=(1./D')*dr*y/(2*tau).*exp(1./(D'*4*tau)*-y.^2);
counts = f*func_array;
end