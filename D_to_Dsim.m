function Dsim = D_to_Dsim(D)

sigma = 0.023;
tau = 0.021742;
Dbound = D;
Dbound = D + .161.*D.^1.53;
Dsim = 3/2*(Dbound - sigma^2/tau);



end 