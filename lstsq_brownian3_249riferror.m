clc
clear 
close all

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


idx = [2, 2, 2,    2, 2];

% Try these error bounds on this parameter
Derror = 0;
ferror = 2; 
de = 0.008;
minus = true;
minus = false;

limit = .85; %chosen
dr = 0.01; %chosen
edges = (0:dr:limit);
counts_exp = histcounts(exp, edges);
dl = 0.001;


%[.003 .004 .050]
 %[.006 .008]

best = 121.7988;

D=[0.088,0.271,1.346];
f=[0.471,0.486,0.043];


Dfit = D;
ffit = [f(1) f(2)];


onebefore = 2000;
v = 1000;
delta = 0;
%% Optimize +/-1 error bar
for t = 1:100
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
difference = abs(onebefore-v);
if difference < 0.02 && dl == 0.001
    dl = 0.005;   
elseif difference < 0.02 && dl == 0.005
    dl = 0.001;    
end
dl = 0.005;
for i = 1:3
if idx(i) == 1 && i ~= Derror
    Dfit(i) = Dfit(i) - dl;
elseif idx(i) == 2 || i == Derror 
    Dfit(i) = Dfit(i);    
else
    Dfit(i) = Dfit(i) + dl;        
end
end
%^^^^^^^^^^^^^^^^^^^^^^^^
for i = 1:2
if ffit(i) <= 0
    ffit(i) = 0;
elseif idx(i+3) == 1 && i ~= ferror
    ffit(i) = ffit(i) - dl;
elseif idx(i+3) == 2 || i == ferror
    ffit(i) = ffit(i);    
else
    ffit(i) = ffit(i) + dl;        
end
end

Darray = zeros([3,3]);
farray = zeros([2,3]);
for i = 1:3
    if i == Derror && minus==false
        Darray(i,:) = [D(i)+de,D(i)+de,D(i)+de];
    elseif i == Derror && minus==true
        Darray(i,:) = [D(i)-de,D(i)-de,D(i)-de];
    else
        Darray(i,:) = Dfit(i)-dl:dl:Dfit(i)+dl;
    end
end
%^^^^^^^^^^^^^^^^^^^^^^^^
for i = 1:2
    if i == ferror && minus==false
        farray(i,:) = [f(i)+de,f(i)+de,f(i)+de];
    elseif i == ferror && minus==true
        farray(i,:) = [f(i)-de,f(i)-de,f(i)-de];
    else
        farray(i,:) = ffit(i)-dl:dl:ffit(i)+dl; 
    end
end

chisq_array = zeros([3,3,3,  3,3]);
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


for m = 1:3
for l = 1:3
for k = 1:3
for b = 1:3
for a = 1:3

nsim=1; %factor to make simulated data more precise
tries = 1;
f3 = 1-farray(1,a)-farray(2,b);


    clc
    disp(t)
    disp(delta)
    fprintf('Fitting distribution for D=[%.3f,%.3f,%.3f], f=[%.3f,%.3f,%.3f]\n\n', ...
        Darray(1,k),Darray(2,l),Darray(3,m),...
        farray(1,a),farray(2,b),f3);

Dsim=[Darray(1,k),Darray(2,l),Darray(3,m)];
fsim=[farray(1,a),farray(2,b),f3];
counts_best = counts_model(Dsim,fsim,limit)*44443;
chisq_array(k,l,m,a,b) = chi_squared(counts_exp, counts_best);



end
end
end
end
end

fprintf('\n');

onebefore = v;

[v, linIdx] = min(chisq_array(:));
[idxC{1:ndims(chisq_array)}] = ind2sub(size(chisq_array),linIdx);
idx = cell2mat(idxC);
delta = v-best; 
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