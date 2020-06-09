function [rne,tt,params,rne_state] = random_diffusion_3D_SphCyl(params)
% To simulate 3D dynamics of RNaseE in 249 and 249-rif
% This is a simplified version of the code used in the paper
% Surovtsev et al. PNAS 2016 113 E7266 
% http://www.pnas.org/content/113/46/E7268.full
%
% Simulations organized as a nested loops with "silent" internal loop that run simulations without saving coordinates, 
%   and external loop in which corrdinates are saved at each iteration.
% All simulations are within a 3D spherocylinder with reflective boundaries  
%
% INPUT:
%   params - parameters for simulation with multiple fields, for example:
%     params.dt - Simulation time step
%     params.dt_out - data output time step
%     params.t_fin - end time of simulation
%   Full list of parameters can be found in a nested function "change_parameters" below
%   Defaults values will be used for parameters value/fields not provided in input
%
% OUTPUT:
%   rne - coordinates of RNaseE vs time, i-th row is a snap shot of coordinates of all RNaseE (X1,Y1, Z1, X2,Y2, Z2 ...) at a given time,(i-1)*dt_out 
%   tt - simulation time
%   params - parameters used for simulation, also contain versio/date info
%   rne_state - state of the RNaseE, i-th row is a list of states (S1,S2,S3,...) at a given time,(i-1)*dt_out, using following convention:
%           0=RNA-bound RNaseE
%           1=free RNaseE
%   
% Andrew Maytin
% 6.8.2020
%__________________________________________________________
%%


%% Parameters
 % First, define all required parameters of simulations 
 % If values provided in the input "params" - change default values to ones provided in params
 % Note: parameters units should be self-consistent. For example, use seconds and micrometers as all time and length units  
 script=mfilename; % save code name for the output 
 params = change_parameters(params); % get all parameters values

% Split 'params' into individual parameters
 % simulation time parameters
dt=params.dt;  % Simulation time step
dt_out=params.dt_out; % data output time step
t_fin=params.t_fin; % end time of simulation
 % simulation box parameters
l0=params.l0;    % cell length
w0=params.w0;  % cell width

 % free and bound RNaseE-related parameters
k_fb=params.k_fb;
k_bf=params.k_bf;
totR=params.totR;  % number of RNaseE
D_bound=params.D_bound; % Diffusion coefficient for bound RNaseE
D_free=params.D_free; % Diffusion coefficient for free RNaseE

%% INTIALIZATION
rng('shuffle'); % to randomly reset RND-generator and change generator method from default

% geometric parameters
l00=l0/2; % long-axis
w00=w0/2;% short axis 
Rcyl=w00; % Radius of the cyliner and caps
Lcyl0=l00-Rcyl; % half-length of cylinder part

% setting 2 states of RNaseE
indxR=1:totR; % just an indexing array 1,2,3,...
boundR=indxR; % initially all ParA are DNA-bound
freeR=[];    % no freely diffusing RNaseE at t=0

% characteristic reaction times, to use in generation of random numbers in stochastic reactions
tau_fb=1/k_fb;  
tau_bf=1/k_bf; 
 
% setting initial coordinates for free/bound RNaseE, randomly distributed   
rand_XYZ=get_rand_XYZ(totR); % get random X, Y and Z-s
x_R0=rand_XYZ(1,:);  
y_R0=rand_XYZ(2,:); 
z_R0=rand_XYZ(3,:); 

x_R=x_R0; 
y_R=y_R0; 
z_R=z_R0; 
 
% making data collectors for the output
rne=111*ones(ceil(t_fin/dt_out)+1,3*totR);         % RNaseE coordinates
rne_state=222*ones(ceil(t_fin/dt_out)+1,totR);    % state of RNaseE
% saving initial values - RNaseE coordinates states - at t=0;  
rne(1,1:3:end)=x_R;
rne(1,2:3:end)=y_R;
rne(1,3:3:end)=z_R;

rne_state(1,:)=zeros(1,totR); % all RNaseE are bound to RNA


%% SIMULATION

ii0=1; % counter for data saving

% initializing list of stochactic times for considered reactions(transitions)
T_fb=[]; % list of stochactis free --> bound transition
T_bf=[]; % list of stochactis bound --> free transition

%pause on

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
for tt1=dt_out:dt_out:t_fin % "save cycle", save data at each step
  
  % Get pools of random times for reaction events
  % if k=0 take a time longer than the total simulation time
  % otherwise draw from exponential distribution with appropriate time constant
  %
  % for free RNaseE -> bound RNaseE transition   
  switch sign(k_fb)
    case 0
      Tpool_fb=10*t_fin*ones([1,2*totR]);  
    case 1
      Tpool_fb=exprnd(tau_fb,[1,2*totR*ceil(dt_out/tau_fb)]);  
  end
  % for bound RNaseE -> free RNaseE transition   
  switch sign(k_bf)
    case 0
      Tpool_bf=10*t_fin*ones([1,2*totR]);  
    case 1
      Tpool_bf=exprnd(tau_bf,[1,2*totR*ceil(dt_out/tau_bf)]);  
  end

 
  % Get pools of random steps, i.e., Brownian dynamics steps
  % Used the formula from https://www.nature.com/articles/nmeth.2367
  Dpool_bound=sqrt(2*D_bound*dt)*normrnd(0,1,[1,3*totR*ceil(dt_out/dt)]); 
  Dpool_free=sqrt(2*D_free*dt)*normrnd(0,1,[1,3*totR*ceil(dt_out/dt)]); 
 
  % counters for drawing numbers from the pools 
  ddB=0; % counter for bound RNaseE diffusion steps
  ddF=0; % counter for free RNaseE diffusion steps
  fb=0;  % counter for free-bound transition
  bf=0;  % counter for bound-free transition
  
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for tt2=dt:dt:dt_out % "silent" cycle, move without saving data
      
    tt=tt1+tt2; % current time
        
    %~~~~~ UPDATE current values:
    % - update lists of Reaction Times since dt time past 
    T_fb=T_fb-dt; 
    T_bf=T_bf-dt; 

    % - update states of ParA
    freeR_new=boundR(T_bf<=0);        % make new transitions: bound -> free
    boundR_new=freeR(T_fb<=0);         % make new transitions: free -> bound
    freeR=[freeR(T_fb>0),freeR_new];  % combine "old" with "new" free RNaseE
    boundR=[boundR(T_bf>0),boundR_new];            % combine "old" with "new" bound RNaseE
    
    % - remove/add Reaction Times according changes in RNaseE states

    % free -> bound transition
    T_fb=T_fb(T_fb>0) ;
    d_indx=length(boundR_new);  % this is how many new bound RNaseE were generated
     if d_indx>0   % if there are some new bound RNaseE
       T_fb=[T_fb, Tpool_fb(fb+1:fb+d_indx)];  
       fb=fb+d_indx; 
     end  
    % bound -> free
    T_bf=T_bf(T_bf>0) ;
    h_indx=length(freeR_new);  % this is how many new free RNaseE were generated
     if h_indx>0  % if there are some new free RNaseE
       T_bf=[T_bf, Tpool_bf(bf+1:bf+h_indx)];  
       bf=bf+h_indx; 
     end  
    %~~~~~ END of "Update current values" 
          
    % ~~~~~ Make BD moves    
 
    % - bound RNaseE moves
       xB_new=x_R(boundR)+Dpool_bound(1,ddB+1:ddB+d_indx);
       yB_new=y_R(boundR)+Dpool_bound(1,ddB+d_indx+1:ddB+2*d_indx);
       zB_new=z_R(boundR)+Dpool_bound(1,ddB+2*d_indx+1:ddB+3*d_indx);
       % apply reflecting boundaries, if necessary,
       [x_R(boundR),y_R(boundR),z_R(boundR)] = apply_boundaries(xB_new,yB_new,zB_new);  
     
     % shift the counter 
     ddB=ddB+3*d_indx;
     
     % - free RNaseE moves
       xF_new=x_R(freeR)+Dpool_free(1,ddF+1:ddF+d_indx);
       yF_new=y_R(freeR)+Dpool_free(1,ddF+d_indx+1:ddF+2*d_indx);
       zF_new=z_R(freeR)+Dpool_free(1,ddF+2*d_indx+1:ddF+3*d_indx);
       % apply reflecting boundaries, if necessary,
       [x_R(freeR),y_R(freeR),z_R(freeR)] = apply_boundaries(xF_new,yF_new,zF_new);  

     % shift the counter 
     ddF=ddF+3*d_indx;
           
  end % "silent" cycle
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   % Current Data collection
   ii0=ii0+1;
   rne(ii0,1:3:end)=x_R;
   rne(ii0,2:3:end)=y_R;
   rne(ii0,3:3:end)=z_R;
   rne_state(ii0,boundR)=0*ones(1,length(boundR));     % 0= bound RNaseE
   rne_state(ii0,freeR)=ones(1,length(freeR));         % 1= free RNaseE

end % "save" cycle
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

%% NESTED FUNCTIONS
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function rand_XYZ=get_rand_XYZ(nn)
% returns n random XYZ triplets ([x1,x2,x3...; y1,y2,y3...; z1,z2,z3...]) within spherocylinder
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

end

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function params_out = change_parameters(params_in)
% declare default values of parameters and change them to values specified in params_in 
  % declare default parameters values
 params0.dt=0.1;  % Simulation time step
 params0.dt_out=1; % data output time step
 params0.t_fin=10; % end time of simulation

 l_0=2.3; w_0=0.5; % cell length and width
 params0.l0=l_0; %l_0=2.5;    % cell length
 params0.w0=w_0; %w_0=0.5;  % cell width

 params0.totR=90;  % number of RNaseE
 %params0.f_bound; % occupancy of bound state

 params0.D_free=0.01; % D coefficient for free RNaseE
 params0.D_bound=0.001; % D coefficient for RNA-bound RNaseE

 params0.k_fb=0.1; % rate of free->bound
 params0.k_bf=0.1; % rate of bound->free
  
 params0.dim=3; % Dimension of simulations;
 params0.script=script; %'diffusing_RNaseE';
 params0.model='2-state model';
 params0.version='2020.06.08';
  
 % check what is provided as input parameters, i.e., in 'params_in', and substitute those in default values
 fldnms_in = fieldnames(params_in);
 fldnms_0 = fieldnames(params0);
 params_out=params0; 
 for kk = 1:length(fldnms_in)
    not_found=1;
    for jj = 1:length(fldnms_0)
      if (strcmpi(fldnms_in{kk},fldnms_0{jj}))
         params_out.(fldnms_0{jj}) = params_in.(fldnms_in{kk});
         not_found=0;
      end
    end
    if not_found
        disp(['Warning: parameter ',fldnms_in{kk},' not found in the list of parameters, check spelling...'])
    end
 end
 if params_out.dt_out<params_out.dt
    params_out.dt_out=params_out.dt; 
 end
  
end


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function [x_new,y_new,z_new] = apply_boundaries(x,y,z)  
% applies reflective boundaries shaped as spherocylinder with cylinder length Lcyl and cylinder and caps radious Rcyl
   x_new=x;
   y_new=y;
   z_new=z;

  xL=abs(x)-Lcyl0;
  %r_old=sqrt(heaviside(xL).*xL.^2+y.^2+z.^2);
  r_old=sqrt((xL>0).*xL.^2+y.^2+z.^2); % using logical seems faster than heaviside
  ind=(r_old>Rcyl);
  r_new=2*Rcyl-r_old(ind);
  x_new(ind)=(xL(ind)>0).*sign(x(ind)).*(Lcyl0+xL(ind).*r_new./r_old(ind)) + (xL(ind)<=0).*x(ind);
  y_new(ind)=y(ind).*r_new./r_old(ind);
  z_new(ind)=z(ind).*r_new./r_old(ind);
end

%% END of ENDS

end