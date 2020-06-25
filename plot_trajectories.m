clc
clear all
close all

%Plot generated RNaseE trajectories

%Set parameters
params.dt=0.1;
params.dt_out=1;
params.t_fin=10;
%params.l0=2.3;
%params.totR=100;
params.D=0.2;

%Simulate trajectories
%    rne - coordinates of RNaseE vs time, i-th row is a snap shot of coordinates of all 
%          RNaseE (X1,Y1, Z1, X2,Y2, Z2 ...) at a given time,(i-1)*dt_out
[rne,tt,params]=rand_diff_3D_SphCyl_theo(params);
x=rne(:,1:3:end);
y=rne(:,2:3:end);
z=rne(:,3:3:end);

%Plot 2D trajectories
plot(x,y)
%hold on
%plot(mean(x),mean(y),'o')

%Plot 3D trajectories
%plot3(x,y,z)

grid on
daspect([1 1 1])

%% Uniform Distribution

dt=params.dt;  % Simulation time step
dt_out=params.dt_out; % data output time step
t_fin=params.t_fin; % end time of simulation
 % simulation box parameters
l0=params.l0;    % cell length
w0=params.w0;  % cell width

 % free and bound RNaseE-related parameters
totR=params.totR;  % number of RNaseE
D=params.D; % Diffusion coefficient

% geometric parameters
l00=l0/2; % long-axis
w00=w0/2;% short axis 
Rcyl=w00; % Radius of the cyliner and caps
Lcyl0=l00-Rcyl; % half-length of cylinder part
%Plot random points in spherocylinder
 nn = 100;
 uni_rand=rand(3,(3+(nn<100)*100)*nn);
 rand_0x=l0*uni_rand(1,:)-l00;
 rand_0y=w0*uni_rand(2,:)-w00;
 rand_0z=w0*uni_rand(3,:)-w00;
 xL=abs(rand_0x)-Lcyl0;
 r_old=sqrt(heaviside(xL).*xL.^2+rand_0y.^2+rand_0z.^2);
 rand_1x=rand_0x(r_old<Rcyl);
 rand_1y=rand_0y(r_old<Rcyl);
 rand_1z=rand_0z(r_old<Rcyl);
 rand_XYZ(1,:)=rand_1x(1:nn);
 rand_XYZ(2,:)=rand_1y(1:nn);
 rand_XYZ(3,:)=rand_1z(1:nn);
 
 hold on
 plot3(rand_1x(1:nn),rand_1y(1:nn),rand_1z(1:nn),'o')
 grid on
 daspect([1 1 1])
