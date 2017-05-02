function phi = BendingModeFun(s,ModeNr)

% Compute clamped-free bending modes for Modenr<5
temp = [1.8751,4.6941,7.8548,10.9955];
lambda= temp(ModeNr);
F = sinh(lambda)+sin(lambda);
G = cosh(lambda)+cos(lambda);
R = s(end);
c1 = G/F ;
phi = cosh(s.*lambda/R)-cos(s.*lambda/R)-c1*(sinh(s.*lambda/R)-sin(s.*lambda/R));