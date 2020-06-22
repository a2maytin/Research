%Averaged positions decreases the variance of single step displacements.
%(by a factor of 1.6)
%
%
%
clc
close all
clear all

% Create simulated trajectories
Dsim = 0.2;
steps = 100;
msteps = 100;
tausim = 0.02;

dt=tausim/msteps;  
dt_out=tausim; %time between each frame (exposure time) in s
t_fin=tausim*steps; 
totR=100;
D=Dsim; %in um^s/s

dT=dt_out;

x_R0=zeros(1,totR);
y_R0=zeros(1,totR);
z_R0=zeros(1,totR);

x_R=x_R0; 
y_R=y_R0; 
z_R=z_R0; 
 
rne=111*ones(floor(t_fin/dt_out)+1,3*totR);  
rne2=111*ones(floor(t_fin/dt_out)+1,3*totR);     

rne(1,1:3:end)=x_R;
rne(1,2:3:end)=y_R;
rne(1,3:3:end)=z_R;

rne2(1,1:3:end)=x_R;
rne2(1,2:3:end)=y_R;
rne2(1,3:3:end)=z_R;

ii0=1; % counter for data saving


% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
for tt1=dt_out:dt_out:t_fin % "save cycle", save data at each step
 
  % Get pool of random steps, i.e., Brownian dynamics steps
  % Used the formula from https://www.nature.com/articles/nmeth.2367
  Dpool=sqrt(2*D*dt)*normrnd(0,1,[1,3*totR*ceil(dt_out/dt)]);  

  % counter for drawing numbers from the pool
  ddB=0; 
  
  % arrays to store the microtrajectories
  x_M = ones(dt_out/dt,totR);
  y_M = ones(dt_out/dt,totR);
  z_M = ones(dt_out/dt,totR);
  
  % counter for storing microtrajectories
  ii1 = 0;
  
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for tt2=dt:dt:dt_out % "silent" cycle, move without saving data
      
    tt=tt1+tt2; % current time
          
    % ~~~~~ Make BD moves    
 
    % - RNaseE moves
       xB_new=x_R+Dpool(1,ddB+1:ddB+totR);
       yB_new=y_R+Dpool(1,ddB+totR+1:ddB+2*totR);
       zB_new=z_R+Dpool(1,ddB+2*totR+1:ddB+3*totR);
       % apply reflecting boundaries
       %[x_R,y_R,z_R] = apply_boundaries(xB_new,yB_new,zB_new);  
       x_R=xB_new;
       y_R=yB_new;
       z_R=zB_new;
       
       ii1 = ii1+1;
       x_M(ii1,:) = x_R;
       y_M(ii1,:) = y_R;
       z_M(ii1,:) = z_R;
     
     % shift the pool counter 
     ddB=ddB+3*totR;
           
  end % "silent" cycle
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   % Current Data collection
   ii0=ii0+1;
   
   rne(ii0,1:3:end)=mean(x_R,1);
   rne(ii0,2:3:end)=mean(y_R,1);
   rne(ii0,3:3:end)=mean(z_R,1);
   
   rne2(ii0,1:3:end)=mean(x_M,1); %actually taking the mean
   rne2(ii0,2:3:end)=mean(y_M,1);
   rne2(ii0,3:3:end)=mean(z_M,1);

end % "save" cycle



h = 1; %counter for pooling displacements

r_model=ones(1,(size(rne,1)-1)*totR);  %array will store simulated displacements

for j = 1:totR %loop through every track
    for i = 1:(size(rne,1)-1) %loop through every x displacement
        r_new = rne(i+1,3*j-2)-rne(i,3*j-2);
        r_model(h) = r_new;
        r_new2 = rne2(i+1,3*j-2)-rne2(i,3*j-2);
        r_model2(h) = r_new2;
        h=h+1;
    end
end

figure
h2 = histogram(r_model)
h2.Normalization = 'probability';
h2.BinWidth = 0.01;
hold on
h = histogram(r_model2)
h.Normalization = 'probability';
h.BinWidth = 0.01;
%just the displacements between t=0 and t=dt_out from averaging (rne2)
r_model3=r_model2(1:(size(rne,1)-1):length(r_model2));
h3 = histogram(r_model3)
h3.Normalization = 'probability';
h3.BinWidth = 0.01;

%plot a gaussian with stdev sqrt(2Dtau)
hold on
y = -0.4:0.01:0.4;
k = sqrt(2 * Dsim * tausim);
k2 = sqrt(2 * Dsim * tausim/1.5);
k3 = sqrt(2 * Dsim * tausim/3);
func = 0.01*normpdf(y,0,k);
plot(y,func,'LineWidth',1.5)
func2 = 0.01*normpdf(y,0,k2);
plot(y,func2,'LineWidth',1.5)
func3 = 0.01*normpdf(y,0,k3);
plot(y,func3,'LineWidth',1.5)
title('Distribution of 1D displacements')
legend('no averaging','averaging','averaging, between t=0 and t=dt')

%plot the tracks
SPACE_UNITS = 'µm';
TIME_UNITS = 's';
tracks = cell(totR, 1);

figure
for i = 1 : totR 

    % Time
    time = (0 : size(rne,1)-1)' * dt_out;

    % Position
    X = [rne(:,3*i-2)];
    X2 = [rne2(:,3*i-2)];

    % Plot
    hold on
    plot(time,X,'b')
    plot(time,X2,'r')

end
clear i X time

