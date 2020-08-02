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


best = 107.9727-0.0243 -0.0747-0.0072 -0.0611 -0.0449-0.0048;
D=[0.071,0.158,0.250,0.469,1.646]
f=[0.283,0.322,0.257,0.110,0.028]




idx = 2*[1, 1, 1, 1, 1,   1, 1, 1,1];

%% Try these error bounds on this parameter
Derror = 0;
ferror = 4; 
de = 0.04;
limit = .85;
dr = 0.01;
edges = (0:dr:limit);
counts_exp = histcounts(exp, edges);
dl = 0.001;

Dfit = D;
ffit = f;
onebefore = 2000;
v = 1000;
delta = 0;
minus = true;
%minus = false;
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

for i = 1:5
if idx(i) == 1 && i ~= Derror
    Dfit(i) = Dfit(i) - dl;
elseif idx(i) == 2 || i == Derror 
    Dfit(i) = Dfit(i);    
else
    Dfit(i) = Dfit(i) + dl;        
end
end
%^^^^^^^^^^^^^^^^^^^^^^^^
for i = 1:4
if ffit(i) <= 0
    ffit(i) = 0;
elseif idx(i+5) == 1 && i ~= ferror
    ffit(i) = ffit(i) - dl;
elseif idx(i+5) == 2 || i == ferror
    ffit(i) = ffit(i);    
else
    ffit(i) = ffit(i) + dl;        
end
end

Darray = zeros([5,3]);
farray = zeros([4,3]);
for i = 1:5
    if i == Derror && minus==false
        Darray(i,:) = [D(i)+de,D(i)+de,D(i)+de];
    elseif i == Derror && minus==true
        Darray(i,:) = [D(i)-de,D(i)-de,D(i)-de];
    else
        Darray(i,:) = Dfit(i)-dl:dl:Dfit(i)+dl;
    end
end
%^^^^^^^^^^^^^^^^^^^^^^^^
for i = 1:4
    if i == ferror && minus==false
        farray(i,:) = [f(i)+de,f(i)+de,f(i)+de];
    elseif i == ferror && minus==true
        farray(i,:) = [f(i)-de,f(i)-de,f(i)-de];
    else
        farray(i,:) = ffit(i)-dl:dl:ffit(i)+dl; 
    end
end

chisq_array = zeros([3,3,3,3,3,  3,3,3,3]);
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


for o = 1:3
for n = 1:3
for m = 1:3
for l = 1:3
for k = 1:3
for d = 1:3
for c = 1:3
for b = 1:3
for a = 1:3

nsim=1; %factor to make simulated data more precise
tries = 1;
f5 = 1-farray(1,a)-farray(2,b)-farray(3,c)-farray(4,d);

if f5 <=0 || Darray(1,k)>Darray(2,l) || Darray(2,l) > Darray(3,m) ...
        || Darray(3,m) > Darray(4,n) || Darray(4,n) > Darray(5,o)
    clc
    disp(t)
    disp(delta)
    chisq_array(k,l,m,n,o,a,b,c,d) = Inf;
    
else
clc
disp(t)
disp(delta)
fprintf('Fitting distribution for D=[%.3f,%.3f,%.3f,%.3f,%.3f], f=[%.3f,%.3f,%.3f,%.3f,%.3f]\n\n', ...
Darray(1,k),Darray(2,l),Darray(3,m),Darray(4,n),Darray(5,o),...
farray(1,a),farray(2,b),farray(3,c),farray(4,d),f5);

Dsim=[Darray(1,k),Darray(2,l),Darray(3,m),Darray(4,n),Darray(5,o)];
fsim=[farray(1,a),farray(2,b),farray(3,c),farray(4,d),f5];
counts_best = counts_model(Dsim,fsim,limit)*44443;
chisq_array(k,l,m,n,o,a,b,c,d) = chi_squared(counts_exp, counts_best);

end



end
end
end
end
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

if f5b <= 0
  fprintf('f5 zeroed out.');
  break  
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
end

f5 = 1- farray(1,idx(6))- farray(2,idx(7))- farray(3,idx(8))- farray(4,idx(9));
fprintf('Best was D=[%.3f,%.3f,%.3f,%.3f,%.3f], f=[%.3f,%.3f,%.3f,%.3f,%.3f]\n\n', ...
        Darray(1,idx(1)),Darray(2,idx(2)),Darray(3,idx(3)),Darray(4,idx(4)),Darray(5,idx(5)),...
        farray(1,idx(6)),farray(2,idx(7)),farray(3,idx(8)),farray(4,idx(9)),f5);

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
%chisq = chisq/(length(exp)-7); 
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