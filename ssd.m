function r = ssd(Dsim,simfac,params1)  
% calculates single step displacements from simulated trajectories
steps = 8*simfac; %<<<<<<< I MADE 7*simfac INTO 8*simfac HERE <<<<<<
msteps = 100;
tausim = 0.021742;
sigma = 0.023;
totR = 7*907;

params1.dt=tausim/msteps;  
params1.dt_out=tausim; %time between each frame (exposure time) in s
params1.t_fin=tausim*steps; 
params1.totR = totR;
params1.sigma = sigma;

% if nargin < 3
%     params1.avg = false;
%     params1.loc = false;
%     params1.sigma = sigma;
%     params1.boundary = true;
% end

r=ones(length(Dsim),(steps-1)*totR);  %array will store simulated displacements

for k = 1:length(Dsim)
params1.D=Dsim(k); %in um^s/s
[rne,~]=rand_diff_3D_SphCyl(params1);

h = 1; %counter for pooling displacements

for j = 1:totR %loop through every track
    for i = 2:(size(rne,1)-1) %loop through every displacement    %<<<<<< I MADE 1 INTO 2 HERE <<<<
        r_new = sqrt((rne(i+1,3*j-2)-rne(i,3*j-2))^2+(rne(i+1,3*j-1)-rne(i,3*j-1))^2);
        r(k,h) = r_new;
        h=h+1;
    end
end

end