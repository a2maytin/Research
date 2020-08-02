clc
clear 
close all

% Import the measured trajectories
file = load('SK249_tracksFinal.mat');
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


idx = [2, 2,      2];

%% Try these error bounds on this parameter
Derror = 0;
ferror = 1; 
de = 0.02;

limit = .5; %chosen
dr = 0.01; %chosen
edges = (0:dr:limit);
counts_exp = histcounts(exp, edges);
dl = 0.005;

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%BEST: 2.4836/116.7301
Dbest = [0.112, 0.437];  %[.004,.015] 
fbest = [0.688, .312]; %[.02 ]

Dfit = Dbest;
ffit = [fbest(1)];


onebefore = 2000;
v = 1000;
delta = 0;
%% Optimize +1 error bar
for t = 1:100
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
difference = abs(onebefore-v);
if difference < 0.02 && dl == 0.001
    dl = 0.005;   
elseif difference < 0.02 && dl == 0.005
    dl = 0.001;    
end
dl = 0.005;
for i = 1:2
if idx(i) == 1 && i ~= Derror
    Dfit(i) = Dfit(i) - dl;
elseif idx(i) == 2 || i == Derror 
    Dfit(i) = Dfit(i);    
else
    Dfit(i) = Dfit(i) + dl;        
end
end
%^^^^^^^^^^^^^^^^^^^^^^^^
for i = 1:1
if ffit(i) <= 0
    ffit(i) = 0;
elseif idx(i+2) == 1 && i ~= ferror
    ffit(i) = ffit(i) - dl;
elseif idx(i+2) == 2 || i == ferror
    ffit(i) = ffit(i);    
else
    ffit(i) = ffit(i) + dl;        
end
end

Darray = zeros([2,3]);
farray = zeros([1,3]);
for i = 1:2
    if i == Derror
        Darray(i,:) = [Dbest(i)+de,Dbest(i)+de,Dbest(i)+de];
    else
        Darray(i,:) = Dfit(i)-dl:dl:Dfit(i)+dl;
    end
end
%^^^^^^^^^^^^^^^^^^^^^^^^
for i = 1:1
    if i == ferror
        farray(i,:) = [fbest(i)+de,fbest(i)+de,fbest(i)+de]; 
    elseif ffit(i) == 0
        farray(i,:) = [0,dl,dl];
    else
        farray(i,:) = ffit(i)-dl:dl:ffit(i)+dl; 
    end
end

chisq_array = zeros([3,3,  3]);
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


for l = 1:3
for k = 1:3
for a = 1:3

nsim=1; %factor to make simulated data more precise
tries = 1;
f2 = 1-farray(1,a);


    clc
    disp(t)
    disp(delta)
    fprintf('Fitting distribution for D=[%.3f,%.3f], f=[%.3f,%.3f]\n\n', ...
        Darray(1,k),Darray(2,l),...
        farray(1,a),f2);

Dsim=[Darray(1,k),Darray(2,l)];
fsim=[farray(1,a),f2];
counts_best = counts_model(Dsim,fsim,limit)*17974;
chisq_array(k,l,a) = chi_squared(counts_exp, counts_best);



end
end
end

fprintf('\n');

onebefore = v;

[v, linIdx] = min(chisq_array(:));
[idxC{1:ndims(chisq_array)}] = ind2sub(size(chisq_array),linIdx);
idx = cell2mat(idxC);
delta = v-116.7301; 
disp(delta)
disp(idx)

if (nnz(2-idx) == 1) || difference <= 1e-3
    fprintf('\nMinimum Found.'); 
    break
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
end
%% Optimize -1 error bar
for t = 1:100
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
difference = abs(onebefore-v);
if difference < 0.02 && dl == 0.001
    dl = 0.005;   
elseif difference < 0.02 && dl == 0.005
    dl = 0.001;    
end
dl = 0.005;
for i = 1:2
if idx(i) == 1 && i ~= Derror
    Dfit(i) = Dfit(i) - dl;
elseif idx(i) == 2 || i == Derror 
    Dfit(i) = Dfit(i);    
else
    Dfit(i) = Dfit(i) + dl;        
end
end
%^^^^^^^^^^^^^^^^^^^^^^^^
for i = 1:1
if ffit(i) <= 0
    ffit(i) = 0;
elseif idx(i+2) == 1 && i ~= ferror
    ffit(i) = ffit(i) - dl;
elseif idx(i+2) == 2 || i == ferror
    ffit(i) = ffit(i);    
else
    ffit(i) = ffit(i) + dl;        
end
end

Darray = zeros([2,3]);
farray = zeros([1,3]);
for i = 1:2
    if i == Derror
        Darray(i,:) = [Dbest(i)-de,Dbest(i)-de,Dbest(i)-de];
    else
        Darray(i,:) = Dfit(i)-dl:dl:Dfit(i)+dl;
    end
end
%^^^^^^^^^^^^^^^^^^^^^^^^
for i = 1:1
    if i == ferror
        farray(i,:) = [fbest(i)-de,fbest(i)-de,fbest(i)-de]; 
    elseif ffit(i) == 0
        farray(i,:) = [0,dl,dl];
    else
        farray(i,:) = ffit(i)-dl:dl:ffit(i)+dl; 
    end
end

chisq_array = zeros([3,3,  3]);
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


for l = 1:3
for k = 1:3
for a = 1:3

nsim=1; %factor to make simulated data more precise
tries = 1;
f2 = 1-farray(1,a);


    clc
    disp(t)
    disp(delta)
    fprintf('Fitting distribution for D=[%.3f,%.3f], f=[%.3f,%.3f]\n\n', ...
        Darray(1,k),Darray(2,l),...
        farray(1,a),f2);

Dsim=[Darray(1,k),Darray(2,l)];
fsim=[farray(1,a),f2];
counts_best = counts_model(Dsim,fsim,limit)*17974;
chisq_array(k,l,a) = chi_squared(counts_exp, counts_best);



end
end
end

fprintf('\n');

onebefore = v;

[v, linIdx] = min(chisq_array(:));
[idxC{1:ndims(chisq_array)}] = ind2sub(size(chisq_array),linIdx);
idx = cell2mat(idxC);
delta = v-116.7301; 
disp(delta)
disp(idx)

if (nnz(2-idx) == 1) || difference <= 1e-3
    fprintf('\nMinimum Found.'); 
    break
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
end
%%

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
%chisq = chisq/(length(exp)-9); 
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function counts = counts_model(D, f, limit)
tau = 0.021742;
%edges = [0:0.016:0.8];
%y = [0.008:0.016:0.792];
dr = 0.01; %chosen
y = [dr/2:dr:limit-dr/2];
func_array=(1./D')*dr*y/(2*tau).*exp(1./(D'*4*tau)*-y.^2);
counts = f*func_array;
end