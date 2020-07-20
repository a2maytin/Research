function Dsim = D_to_Dsim(D)

sigma = 0.04;
tau = 0.021742;
Dbound = D;
Dbound = D + D.^2.*0.19+0.0;

Dsim = 3/2*(Dbound - sigma^2/tau);



end