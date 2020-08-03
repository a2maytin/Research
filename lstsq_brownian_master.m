% Perform least squares fit to experimental data.
% Find optimal values of the fit parameters D and f.

% This code iteratively finds the value of the parameters that minimizes 
% the sum of square errors. Note: Usually least squares fitting is done 
% with the Levenberg Marquadt Algorithm which uses gradient descent. This 
% code does not use gradient descent, instead scanning all possible parameter
% combinations step size = dl away from the current set of parameters.

% Fit judged poor when chi2 is greater than best fit chi2 by one degree of 
% freedom. Therefore, the uncertainty in a fit parameter (f) is taken to be 
% the deviation (de) in that fit parameter such when that parameter is 
% fixed at p±de and all other parameters are optimized, results in a chi2 
% of (bestfit chi2) +1. Source:
% https://www.phys.hawaii.edu/~varner/PHYS305-Spr12/DataFitting.html

% Import the measured trajectories
file = load('SK249_tracksFinal.mat');  %SK249-rif.tracksFinal.mat
traj = file.tracksFinal;
coord = {traj.tracksCoordAmpCG};
pos = {traj.tracksCoordXY};

% convert units to SI units, specific to 249/249-rif data for this project
pixelSize = .16; % in micrometers
timeStep = 21.742e-3;

exp=ones(1,(sum(cellfun(@length, pos))-numel(pos))); %store displacements

k = 1; %counter for pooling displacements
for i = 1:numel(pos) %loop through every track
    for j = 1:(length(pos{i})-1) %loop through every displacement
        r_new = sqrt((pos{i}(j+1,1)-pos{i}(j,1))^2+(pos{i}(j+1,2)-pos{i}(j,2))^2);
        exp(k) = r_new*pixelSize;
        k=k+1;
    end
end
%% INITIAL GUESS, SET FIXED PARAMETERS

numStates = 5; % Number of states, can be 1-5

limit = .8; % Set the limit for binning to avoid zero count bins
dr = 0.01; % Set the bin width
dl = 0.001; % parameter increment used during optimization
maxSteps = 100; % Maximum number of iterations before optimization times out 

% Initial Guess
D0=[.038, .1, .19, .56, 1.6]; 
f0=[.074, .28, .48, .15]; 

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% dirStep : Descent direction
% 1 = step down
% 2 = stay same
% 3 = step up
dirStep = 2*ones([1,9]); 
edges = (0:dr:limit);
counts_exp = histcounts(exp, edges);
Dfit = D0;
ffit = f0;
N = length(exp); 
chisq_history = Inf*ones([1,maxSteps+1]);

%% ITERATION: FIND BEST FIT PARAMETERS
for t = 1:maxSteps
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% UPDATE D AND f
for i = 1:numStates
if dirStep(i) == 1 
    Dfit(i) = Dfit(i) - dl;
elseif dirStep(i) == 2
    Dfit(i) = Dfit(i);    
else
    Dfit(i) = Dfit(i) + dl;        
end
end
for i = 1:numStates-1
if ffit(i) <= 0 % restrict f = [0,1)
    ffit(i) = 0;
elseif dirStep(i+5) == 1 
    ffit(i) = ffit(i) - dl;
elseif dirStep(i+5) == 2 
    ffit(i) = ffit(i);    
else
    ffit(i) = ffit(i) + dl;        
end
end
% Parameter values to be tested during this iteration
Darray = Inf*ones([5,3]);
farray = zeros([5,3]);
for i = 1:numStates
        Darray(i,:) = [Dfit(i)-dl Dfit(i) Dfit(i)+dl];
end
for i = 1:numStates-1
        farray(i,:) = [ffit(i)-dl ffit(i) ffit(i)+dl]; 
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% MAIN LOOP 
chisq_array = Inf*ones(3*ones([1,9])); % Array to store the sum of square resdiuals (chi-squared value)
% Loop through D
for o = 1:(1+2*(numStates>=5))
for n = 1:(1+2*(numStates>=4))
for m = 1:(1+2*(numStates>=3))
for l = 1:(1+2*(numStates>=2))
for k = 1:(1+2*(numStates>=1))
% Loop through f   
for d = 1:(1+2*(numStates>4))
for c = 1:(1+2*(numStates>3))
for b = 1:(1+2*(numStates>2))
for a = 1:(1+2*(numStates>1))

Dind = [k,l,m,n,o];
find = [a,b,c,d];
farray(numStates,1) = 1;
for r = 1:numStates-1
farray(numStates,1) = farray(numStates,1) - farray(r,find(r)); % Get the last occupancy
end

if farray(numStates,1) <0 || Darray(1,k)>Darray(2,l) || Darray(2,l) > Darray(3,m) ...
        || Darray(3,m) > Darray(4,n) || Darray(4,n) > Darray(5,o)
    chisq_array(k,l,m,n,o,a,b,c,d) = Inf; %Invalid Fit, set chisq = Inf
else

Dsim=[Darray(1,k),Darray(2,l),Darray(3,m),Darray(4,n),Darray(5,o)];
fsim=[farray(1,a),farray(2,b),farray(3,c),farray(4,d),farray(5,1)];
counts_best = counts_model(Dsim,fsim,timeStep, dr, limit)*N;
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
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

% FIND MINIMUM THIS STEP
[chisq_min, linIdx] = min(chisq_array(:));
[idxC{1:ndims(chisq_array)}] = ind2sub(size(chisq_array),linIdx);
dirStep = cell2mat(idxC);

chisq_history(t+1) = chisq_min;
improvement = chisq_history(t) - chisq_history(t+1);
farray(numStates,1) = 0;
fN = 1- farray(1,dirStep(6))- farray(2,dirStep(7))- farray(3,dirStep(8))- farray(4,dirStep(9));

% Print smallest chisquared this step
Dbest = Darray(1,dirStep(1));
for x = 2:numStates
Dbest = [Dbest , Darray(x,dirStep(x))];    
end
fbest = farray(1,dirStep(1));
for x = 2:numStates-1
fbest = [fbest , farray(x,dirStep(x))];    
end
fbest = [fbest, fN];
msg = strcat('Lowest Chi-squared this Step=', num2str(chisq_min) ,  ', found at:\n D = [' , num2str(Dbest) , ']\n f = [' , num2str(fbest) , ']\n\n');
clc
fprintf(msg)
if (nnz(2-dirStep) == 0) || improvement <= 1e-3
    fprintf('\nMinimum Reached.\n'); 
    chisq_best = chisq_min;
    break
end

end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%% INITIAL GUESS PARAMETER ERROR
% To find the error (uncertainty) in parameter p, the smallest variation de 
% is found such that when p is fixed at p±de and all other parameters are
% optimized, the best fit chi-squared is >= one more than the
% best fit chi-squared. 

% For example, for a two (p1,p2) parameter fit, the grid of chi-squared 
% could looklike this, with best fit chi-squared = 40:

%     \ p2  .01     .02      .03     .04   .05      .06
%  p1  \__________________________________________________
%  .05 |  43.00   41.00   41.00   41.00   41.00   41.00
%  .06 |  41.00   41.00   40.84   40.94   41.00   41.00
%  .07 |  41.00   40.33   40.13   40.23   41.00   41.00
%  .08 |  41.00   40.20  (40.00)  40.10   41.00   41.00
%  .09 |  41.00   40.37   40.17   40.26   41.00   41.00
%  .10 |  41.00   40.87   40.67   40.77   41.00   41.00
%  .11 |  41.00   41.00   41.00   41.00   41.00   43.00  

% Examining the chi-squared = 41 countour, we can report the best fit parameters (p1,p2) = .08±.03, .03±.02

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Derror = 4; %Pick parameter to find error bounds (set to 0 if checking f)
ferror = 0; %Pick parameter to find error bounds (set to 0 if checking D)

de = 0.1; %Try these error bounds on this parameter
minus = true; % True if checking negative error bar and vice versa
maxStepsError = 10;


% Dbest = ; %best fit parameters from above section
% fbest = ;
% chisq_best = ; %best chi-squared from above section
Dfit = Dbest;
ffit = fbest;
minimum_reached = false;
%% FIND ERROR ON BEST FIT PARAMETERS
for t = 1:maxStepsError
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% UPDATE D AND f
for i = 1:numStates
if dirStep(i) == 1 && i ~= Derror
    Dfit(i) = Dfit(i) - dl;
elseif dirStep(i) == 2 == 2 || i == Derror 
    Dfit(i) = Dfit(i);    
else
    Dfit(i) = Dfit(i) + dl;        
end
end
for i = 1:numStates-1
if ffit(i) <= 0 % restrict f = [0,1)
    ffit(i) = 0;
elseif dirStep(i+5)  == 1 && i ~= ferror
    ffit(i) = ffit(i) - dl;
elseif dirStep(i+5) == 2 || i == ferror
    ffit(i) = ffit(i);    
else
    ffit(i) = ffit(i) + dl;        
end
end
% Parameter values to be tested during this iteration
Darray = Inf*ones([5,3]);
farray = zeros([5,3]);
for i = 1:numStates
    if i == Derror && minus==false
        Darray(i,:) = [Dbest(i)+de,Dbest(i)+de,Dbest(i)+de];
    elseif i == Derror && minus==true
        Darray(i,:) = [Dbest(i)-de,Dbest(i)-de,Dbest(i)-de];
    else
        Darray(i,:) = Dfit(i)-dl:dl:Dfit(i)+dl;
    end
end
for i = 1:numStates-1
    if i == ferror && minus==false
        farray(i,:) = [fbest(i)+de,fbest(i)+de,fbest(i)+de];
    elseif i == ferror && minus==true
        farray(i,:) = [fbest(i)-de,fbest(i)-de,fbest(i)-de];
    else
        farray(i,:) = ffit(i)-dl:dl:ffit(i)+dl; 
    end
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% MAIN LOOP 
chisq_array = Inf*ones(3*ones([1,9])); % Array to store the sum of square resdiuals (chi-squared value)
% Loop through D
for o = 1:(1+2*(numStates>=5))
for n = 1:(1+2*(numStates>=4))
for m = 1:(1+2*(numStates>=3))
for l = 1:(1+2*(numStates>=2))
for k = 1:(1+2*(numStates>=1))
% Loop through f   
for d = 1:(1+2*(numStates>4))
for c = 1:(1+2*(numStates>3))
for b = 1:(1+2*(numStates>2))
for a = 1:(1+2*(numStates>1))

Dind = [k,l,m,n,o];
find = [a,b,c,d];
farray(numStates,1) = 1;
for r = 1:numStates-1
farray(numStates,1) = farray(numStates,1) - farray(r,find(r)); % Get the last occupancy
end

if farray(numStates,1) <0 || Darray(1,k)>Darray(2,l) || Darray(2,l) > Darray(3,m) ...
        || Darray(3,m) > Darray(4,n) || Darray(4,n) > Darray(5,o)
    chisq_array(k,l,m,n,o,a,b,c,d) = Inf; %Invalid Fit, set chisq = Inf
else

Dsim=[Darray(1,k),Darray(2,l),Darray(3,m),Darray(4,n),Darray(5,o)];
fsim=[farray(1,a),farray(2,b),farray(3,c),farray(4,d),farray(5,1)];
counts_best = counts_model(Dsim,fsim,timeStep, dr, limit)*N;
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
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

% FIND MINIMUM THIS STEP
[chisq_min, linIdx] = min(chisq_array(:));
[idxC{1:ndims(chisq_array)}] = ind2sub(size(chisq_array),linIdx);
dirStep = cell2mat(idxC);

chisq_history(t+1) = chisq_min;
improvement = chisq_history(t) - chisq_history(t+1);

farray(numStates,1) = 0;
fN = 1- farray(1,dirStep(6))- farray(2,dirStep(7))- farray(3,dirStep(8))- farray(4,dirStep(9));

% Print smallest chisquared this step
DbestErr = Darray(1,dirStep(1));
for x = 2:numStates
DbestErr = [DbestErr , Darray(x,dirStep(x))];    
end
fbestErr = farray(1,dirStep(1));
for x = 2:numStates-1
fbestErr = [fbestErr , farray(x,dirStep(x))];    
end
fbestErr = [fbestErr, fN];
msg = strcat('Lowest Chi-squared this Step=', num2str(chisq_min) ,  ', found at:\n D = [' , num2str(DbestErr) , ']\n f = [' , num2str(fbestErr) , ']\n\n');
clc
fprintf(msg)
if (nnz(2-dirStep) == 0) || improvement <= 1e-4
    fprintf('\nMinimum Reached.\n'); 
    minimum_reached = true;
    break
end

end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
if abs(chisq_min - chisq_best) < .7
    fprintf('The error bar of %.3f is too low, try larger error bar',de)
elseif abs(chisq_min - chisq_best) > .7 && abs(chisq_min - chisq_best) < 1
    fprintf('The error bar of %.3f may be slightly too low,\nbecause it yields a chi-squared %.2f higher than best fit',de,abs(chisq_min - chisq_best))
elseif minimum_reached == true && abs(chisq_min - chisq_best) > 2
    fprintf('The error bar of %.3f is too high, try smaller error bar',de) 
elseif minimum_reached == true && abs(chisq_min - chisq_best) > 1.1 && abs(chisq_min - chisq_best) < 2
    fprintf('The error bar of %.3f may be a slightly high estimate,\nbecause it yields a chi-squared %.2f higher than best fit',de,abs(chisq_min - chisq_best)) 
elseif minimum_reached == true && abs(chisq_min - chisq_best) > 1 && abs(chisq_min - chisq_best) < 1.1
    fprintf('The error bar of %.3f is a great estimate,\nbecause it yields a chi-squared %.2f higher than best fit',de,abs(chisq_min - chisq_best)) 
elseif minimum_reached == false 
    fprintf('Minimum not reached, run section again or try smaller error bar') 
end
fprintf('\n')

%% NATIVE FUNCTIONS
function chisq = chi_squared(exp,sim)  
residuals = exp-sim;
% assume Poisson statistics
errors = max(1,sqrt(exp));
pulls = residuals./errors;
chisq = sum(pulls.*pulls);
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function counts = counts_model(D, f, tau, dr, limit)
y = (dr/2:dr:limit-dr/2);
func_array=(1./D')*dr*y/(2*tau).*exp(1./(D'*4*tau)*-y.^2);
counts = f*func_array;
end