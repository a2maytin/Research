function r = ssd(Dsim,N,params1)  
% Convenience function that calculates single step displacements from 
% simulated trajectories created by rand_diff_3D-SphyCyl, for use in 
% lstsq_simulation

% INPUT:
%   Dsim = simulated diffusion coefficient
%   N - goal number of total displacements
%   params1 - parameters for simulation
%
% OUTPUT:
%   r - pool of single step displacements from the simulated data
%   NOTE: When getting/using data from simulation, I discard position at time 
%   t = 0 (initial position), because it takes into account neither 
%   averaging due to finite camera exposure nor localization error.
%__________________________________________________________________________

steps = 8;          %total useful frames+1
msteps = 100;       %number of microsteps
tausim = 0.021742;  %time step 
sigma = 0.023;      %localization error from UTrack data structure 
totR = ceil(N/(steps-1));

% set the missing fields 
params1.dt=tausim/msteps;  
params1.dt_out=tausim; %time between each frame (exposure time) in s
params1.t_fin=tausim*steps; 
params1.totR = totR;
params1.sigma = sigma;

r=ones(length(Dsim),(steps-1)*totR);  %array will store simulated displacements

for k = 1:length(Dsim)
params1.D=Dsim(k); %in um^s/s
[rne,~]=rand_diff_3D_SphCyl(params1);

h = 1; %counter for pooling displacements

for j = 1:totR %loop through every track
    for i = 2:(size(rne,1)-1) %loop through every displacement (except for first one, which is not useful)
        r_new = sqrt((rne(i+1,3*j-2)-rne(i,3*j-2))^2+(rne(i+1,3*j-1)-rne(i,3*j-1))^2);
        r(k,h) = r_new;
        h=h+1;
    end
end

end