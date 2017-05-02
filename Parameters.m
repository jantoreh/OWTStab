function par = Parameters
% Rotor parameters
par.R    = 50;         % [m]
par.c    = 3;          % [m]
par.acg  = 1.2;        % [m]
par.EIx  = 9.87975e9;  % [Nm^2]
par.EIy  = 17.56404e9; % [Nm^2]
par.GK   = 0.1764804e9;% [Nm^2]
par.m    = 220;        % [kg/m]
par.J    = 275.75;     % [khm^2]
% Nacelle parameters
par.Ls   = 5;          % [m]
par.Gs   = 0.5e9;      % [Nm/rad]
par.Mn   = 205e3;      % [kg]
par.Ix   = 4500e3;     % [kgm^2]
par.Iy   = 1200e3;     % [kgm^2]
par.Iz   = 4500e3;     % [kgm^2]
% Damping in logarithmic decrement
par.etaf = 0.01;       % [-]
par.etae = 0.01;       % [-]
par.etab = 0.04;       % [-]
par.etat = 0.005;      % [-]
par.etad = 0.3;        % [-]
% Air density
par.rho  = 1.225;      % [kg/m^3]
% Nacelle/tower top parameters
D    = 5.00;  % [m]
d    = 4.92;  % [m]
Et   = 211e9; % [N/m^2]
ny   = 0.33;  % [-]
rhot = 7850;  % [kg/m^3]
H    = 70;    % [m]
% Tower stiffness parameters derived from elementary loading cases
% [Dubbel, Handbook of Mech. Eng.]
EIt = Et*pi*(D^4-d^4)/64;
GKt = Et/(2*(1+ny))*pi*(D^4-d^4)/32;
par.Kx  =  12*EIt/H^3; 
par.Ky  =  par.Kx;
par.Gx  =  4*EIt/H;
par.Gy  =  par.Gx;
par.Gz  =  GKt/H;
par.gxy = -6*EIt/H^2;
% Nacelle/tower equivalent mass
mt  = pi*(D^2-d^2)/4*rhot;
h   = linspace(0,H,500);
phi = BendingModeFun(h,1);
par.Mtot = par.Mn+trapz(h,phi.^2*mt)/phi(end)^2;