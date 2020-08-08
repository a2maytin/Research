function Dsim = D_to_Dsim(D)
% Take the apparent D values from the brownian diffusion fit, and 
% back-calculate the best D values to use for the simulation.

sigma = 0.024;
tau = 0.021742;

% REVERSE THE MICROSCOPE MEASUREMENT EFFECT:
% The effects of averaging and localization error on the apparent D are
% well documented in Xavier Michalet's paper on the subject:
% Michalet X. (2010). Mean square displacement analysis of single-particle 
% trajectories with localization error: Brownian motion in an isotropic 
% medium. Physical review. E, Statistical, nonlinear, and soft matter 
% physics, 82(4 Pt 1), 041914.

Dbound = 3/2*(D - sigma^2/tau);

% REVERSE THE BOUNDARY EFFECT:
% The numbers .161 and 1.53 were derived by logarithmically fitting
% the effect of the boundary on multiple simulated D values in a spherocylinder
% with length 3.0 microns, width 1.1 microns (you can also look at my 
% 7/20/2020 presentation in the google drive for the relevant plots).
Dsim = Dbound + .161.*Dbound.^1.53;



end 