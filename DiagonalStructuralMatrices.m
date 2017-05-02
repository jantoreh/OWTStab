function [Mb,Kb,Gb,Mt,Gt,Kt] = DiagonalStructuralMatrices(par,Omega)
% Turbine parameters
R    = par.R   ;   
c    = par.c   ;   
a_cg = par.acg ;   
EIx  = par.EIx ;   
EIy  = par.EIy ;   
GK   = par.GK  ;   
m    = par.m   ;   
J    = par.J   ;   
Ls   = par.Ls  ;   
Gs   = par.Gs  ;
Mn   = par.Mn  ;   
Ix   = par.Ix  ;   
Iy   = par.Iy  ;   
Iz   = par.Iz  ;   
M    = par.Mtot;
Kx   = par.Kx  ;
Ky   = par.Ky  ;
Gx   = par.Gx  ;
Gy   = par.Gy  ;
Gz   = par.Gz  ;
gxy  = par.gxy ;
% Blade data
lambda1 = 1.8751;
c1 = .7340955139;
% Blade mass matrix
Mb = [
-m * R * (0.4e1 * exp((2 * lambda1)) * c1 + 0.1e1 - c1 ^ 2 * exp((4 * lambda1)) + 0.4e1 * c1 ^ 2 * cos(lambda1) * exp(lambda1) + 0.4e1 * cos(lambda1) * exp((3 * lambda1)) + 0.4e1 * c1 ^ 2 * sin(lambda1) * exp(lambda1) + 0.4e1 * sin(lambda1) * exp((3 * lambda1)) + 0.4e1 * c1 ^ 2 * exp((2 * lambda1)) * cos(lambda1) * sin(lambda1) + 0.2e1 * c1 * exp((4 * lambda1)) - exp((4 * lambda1)) - 0.4e1 * c1 ^ 2 * exp((3 * lambda1)) * cos(lambda1) - 0.8e1 * lambda1 * exp((2 * lambda1)) + 0.4e1 * c1 ^ 2 * exp((3 * lambda1)) * sin(lambda1) - 0.4e1 * cos(lambda1) * sin(lambda1) * exp((2 * lambda1)) + 0.8e1 * c1 * sin(lambda1) * exp(lambda1) - 0.8e1 * c1 * cos(lambda1) ^ 2 * exp((2 * lambda1)) + c1 ^ 2 - 0.4e1 * cos(lambda1) * exp(lambda1) + 0.4e1 * sin(lambda1) * exp(lambda1) - 0.8e1 * c1 * exp((3 * lambda1)) * sin(lambda1) + 0.2e1 * c1) / lambda1 * exp(-(2 * lambda1)) / 0.8e1 0 0.2e1 * lambda1 * R * a_cg * m / (-pi ^ 4 + (16 * lambda1 ^ 4)) * (-0.8e1 * pi * lambda1 + 0.2e1 * cos(lambda1) * c1 * pi ^ 2 - exp(-lambda1) * c1 * pi ^ 2 + 0.2e1 * sin(lambda1) * pi ^ 2 - 0.4e1 * (lambda1 ^ 2) * exp(lambda1) + 0.4e1 * exp(-lambda1) * (lambda1 ^ 2) * c1 + 0.8e1 * sin(lambda1) * (lambda1 ^ 2) + pi ^ 2 * exp(lambda1) - c1 * pi ^ 2 * exp(lambda1) + 0.8e1 * cos(lambda1) * (lambda1 ^ 2) * c1 - exp(-lambda1) * pi ^ 2 + 0.4e1 * (lambda1 ^ 2) * c1 * exp(lambda1) + 0.4e1 * exp(-lambda1) * (lambda1 ^ 2))
 0 -m * R * (0.4e1 * exp((2 * lambda1)) * c1 + 0.1e1 - c1 ^ 2 * exp((4 * lambda1)) + 0.4e1 * c1 ^ 2 * cos(lambda1) * exp(lambda1) + 0.4e1 * cos(lambda1) * exp((3 * lambda1)) + 0.4e1 * c1 ^ 2 * sin(lambda1) * exp(lambda1) + 0.4e1 * sin(lambda1) * exp((3 * lambda1)) + 0.4e1 * c1 ^ 2 * exp((2 * lambda1)) * cos(lambda1) * sin(lambda1) + 0.2e1 * c1 * exp((4 * lambda1)) - exp((4 * lambda1)) - 0.4e1 * c1 ^ 2 * exp((3 * lambda1)) * cos(lambda1) - 0.8e1 * lambda1 * exp((2 * lambda1)) + 0.4e1 * c1 ^ 2 * exp((3 * lambda1)) * sin(lambda1) - 0.4e1 * cos(lambda1) * sin(lambda1) * exp((2 * lambda1)) + 0.8e1 * c1 * sin(lambda1) * exp(lambda1) - 0.8e1 * c1 * cos(lambda1) ^ 2 * exp((2 * lambda1)) + c1 ^ 2 - 0.4e1 * cos(lambda1) * exp(lambda1) + 0.4e1 * sin(lambda1) * exp(lambda1) - 0.8e1 * c1 * exp((3 * lambda1)) * sin(lambda1) + 0.2e1 * c1) / lambda1 * exp(-(2 * lambda1)) / 0.8e1 0
  0.2e1 * lambda1 * R * a_cg * m / (-pi ^ 4 + (16 * lambda1 ^ 4)) * (-0.8e1 * pi * lambda1 + 0.2e1 * cos(lambda1) * c1 * pi ^ 2 - exp(-lambda1) * c1 * pi ^ 2 + 0.2e1 * sin(lambda1) * pi ^ 2 - 0.4e1 * (lambda1 ^ 2) * exp(lambda1) + 0.4e1 * exp(-lambda1) * (lambda1 ^ 2) * c1 + 0.8e1 * sin(lambda1) * (lambda1 ^ 2) + pi ^ 2 * exp(lambda1) - c1 * pi ^ 2 * exp(lambda1) + 0.8e1 * cos(lambda1) * (lambda1 ^ 2) * c1 - exp(-lambda1) * pi ^ 2 + 0.4e1 * (lambda1 ^ 2) * c1 * exp(lambda1) + 0.4e1 * exp(-lambda1) * (lambda1 ^ 2)) 0 R * m * a_cg ^ 2 / 0.2e1 + R * J / 0.2e1
  ];
% Blade stiffness matrix
Kb = [
exp(-(2 * lambda1)) * (-(12 * EIx * lambda1 ^ 4) - (12 * EIx * lambda1 ^ 4 * c1 ^ 2) + 0.96e2 * EIx * (lambda1 ^ 5) * exp((2 * lambda1)) - (24 * EIx * lambda1 ^ 4 * c1) - 0.48e2 * EIx * (lambda1 ^ 4) * cos(lambda1) * exp(lambda1) + 0.48e2 * EIx * (lambda1 ^ 4) * sin(lambda1) * exp(lambda1) + 0.48e2 * EIx * (lambda1 ^ 4) * cos(lambda1) * sin(lambda1) * exp((2 * lambda1)) + 0.96e2 * EIx * (lambda1 ^ 4) * c1 * cos(lambda1) ^ 2 * exp((2 * lambda1)) + 0.96e2 * EIx * (lambda1 ^ 4) * c1 * sin(lambda1) * exp(lambda1) - 0.48e2 * EIx * (lambda1 ^ 4) * (c1 ^ 2) * exp((2 * lambda1)) * cos(lambda1) * sin(lambda1) - 0.24e2 * EIx * (lambda1 ^ 4) * c1 * exp((4 * lambda1)) - 0.48e2 * EIx * (lambda1 ^ 4) * exp((2 * lambda1)) * c1 + 0.12e2 * EIx * (lambda1 ^ 4) * (c1 ^ 2) * exp((4 * lambda1)) + 0.48e2 * EIx * (lambda1 ^ 4) * (c1 ^ 2) * sin(lambda1) * exp(lambda1) + 0.12e2 * EIx * (lambda1 ^ 4) * exp((4 * lambda1)) + 0.48e2 * EIx * (lambda1 ^ 4) * cos(lambda1) * exp((3 * lambda1)) + 0.48e2 * EIx * (lambda1 ^ 4) * sin(lambda1) * exp((3 * lambda1)) - 0.48e2 * EIx * (lambda1 ^ 4) * (c1 ^ 2) * exp((3 * lambda1)) * cos(lambda1) + 0.48e2 * EIx * (lambda1 ^ 4) * (c1 ^ 2) * cos(lambda1) * exp(lambda1) - 0.96e2 * EIx * (lambda1 ^ 4) * c1 * exp((3 * lambda1)) * sin(lambda1) + 0.48e2 * EIx * (lambda1 ^ 4) * (c1 ^ 2) * exp((3 * lambda1)) * sin(lambda1) + (6 * m * Omega ^ 2 * R ^ 4 * c1 ^ 2 * lambda1) + (12 * m * Omega ^ 2 * R ^ 4 * c1 * lambda1) + (3 * m * Omega ^ 2 * R ^ 4) + (3 * m * Omega ^ 2 * R ^ 4 * c1 ^ 2) + (6 * m * Omega ^ 2 * R ^ 4 * c1) + (6 * m * Omega ^ 2 * R ^ 4 * lambda1) - 0.3e1 * m * (Omega ^ 2) * (R ^ 4) * exp((4 * lambda1)) + 0.32e2 * m * (Omega ^ 2) * (R ^ 4) * (lambda1 ^ 3) * (c1 ^ 2) * exp((2 * lambda1)) + 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * c1 * sin(lambda1) * exp((3 * lambda1)) * lambda1 - 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * cos(lambda1) * lambda1 * exp(lambda1) - 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * (lambda1 ^ 2) * c1 * exp((2 * lambda1)) + 0.6e1 * m * (Omega ^ 2) * (R ^ 4) * exp((4 * lambda1)) * (c1 ^ 2) * lambda1 + 0.24e2 * m * (Omega ^ 2) * (R ^ 4) * cos(lambda1) ^ 2 * lambda1 * exp((2 * lambda1)) - 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * cos(lambda1) * exp((3 * lambda1)) * lambda1 - 0.12e2 * m * (Omega ^ 2) * (R ^ 4) * exp((4 * lambda1)) * c1 * lambda1 + 0.12e2 * m * (Omega ^ 2) * (R ^ 4) * (c1 ^ 2) * exp((2 * lambda1)) * lambda1 - 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * c1 * cos(lambda1) * sin(lambda1) * lambda1 * exp((2 * lambda1)) + 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * c1 * sin(lambda1) * lambda1 * exp(lambda1) + 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * (c1 ^ 2) * sin(lambda1) * lambda1 * exp(lambda1) - 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * c1 * cos(lambda1) * lambda1 * exp(lambda1) + 0.6e1 * m * (Omega ^ 2) * (R ^ 4) * lambda1 * exp((4 * lambda1)) + 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * exp((3 * lambda1)) * c1 * cos(lambda1) * lambda1 - 0.24e2 * m * (Omega ^ 2) * (R ^ 4) * (c1 ^ 2) * exp((2 * lambda1)) * lambda1 * cos(lambda1) ^ 2 - 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * (c1 ^ 2) * sin(lambda1) * exp((3 * lambda1)) * lambda1 + 0.12e2 * m * (Omega ^ 2) * (R ^ 4) * exp((2 * lambda1)) * c1 - 0.3e1 * m * (Omega ^ 2) * (R ^ 4) * (c1 ^ 2) * exp((4 * lambda1)) + 0.24e2 * m * (Omega ^ 2) * (R ^ 4) * exp(lambda1) * (c1 ^ 2) * cos(lambda1) + 0.24e2 * m * (Omega ^ 2) * (R ^ 4) * cos(lambda1) * exp((3 * lambda1)) + 0.24e2 * m * (Omega ^ 2) * (R ^ 4) * (c1 ^ 2) * sin(lambda1) * exp(lambda1) + 0.24e2 * m * (Omega ^ 2) * (R ^ 4) * sin(lambda1) * exp((3 * lambda1)) + 0.12e2 * m * (Omega ^ 2) * (R ^ 4) * (c1 ^ 2) * exp((2 * lambda1)) * cos(lambda1) * sin(lambda1) + 0.6e1 * m * (Omega ^ 2) * (R ^ 4) * c1 * exp((4 * lambda1)) - 0.24e2 * m * (Omega ^ 2) * (R ^ 4) * (c1 ^ 2) * exp((3 * lambda1)) * cos(lambda1) - 0.12e2 * m * (Omega ^ 2) * (R ^ 4) * lambda1 * exp((2 * lambda1)) + 0.24e2 * m * (Omega ^ 2) * (R ^ 4) * (c1 ^ 2) * exp((3 * lambda1)) * sin(lambda1) - 0.12e2 * m * (Omega ^ 2) * (R ^ 4) * cos(lambda1) * sin(lambda1) * exp((2 * lambda1)) + 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * c1 * sin(lambda1) * exp(lambda1) - 0.24e2 * m * (Omega ^ 2) * (R ^ 4) * c1 * cos(lambda1) ^ 2 * exp((2 * lambda1)) - 0.24e2 * m * (Omega ^ 2) * (R ^ 4) * cos(lambda1) * exp(lambda1) + 0.24e2 * m * (Omega ^ 2) * (R ^ 4) * exp(lambda1) * sin(lambda1) - 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * c1 * exp((3 * lambda1)) * sin(lambda1)) / (R ^ 3) / lambda1 / 0.96e2 0 0
 0 exp(-(2 * lambda1)) * ((6 * m * Omega ^ 2 * R ^ 4 * c1 ^ 2 * lambda1) + (12 * m * Omega ^ 2 * R ^ 4 * c1 * lambda1) + (15 * m * Omega ^ 2 * R ^ 4) + (15 * m * Omega ^ 2 * R ^ 4 * c1 ^ 2) + (30 * m * Omega ^ 2 * R ^ 4 * c1) + (6 * m * Omega ^ 2 * R ^ 4 * lambda1) - 0.15e2 * m * (Omega ^ 2) * (R ^ 4) * exp((4 * lambda1)) - 0.48e2 * EIy * (lambda1 ^ 4) * exp((2 * lambda1)) * c1 + 0.12e2 * EIy * (lambda1 ^ 4) * (c1 ^ 2) * exp((4 * lambda1)) - 0.48e2 * EIy * (lambda1 ^ 4) * cos(lambda1) * exp(lambda1) + 0.48e2 * EIy * (lambda1 ^ 4) * exp(lambda1) * sin(lambda1) - 0.24e2 * EIy * (lambda1 ^ 4) * c1 * exp((4 * lambda1)) + 0.48e2 * EIy * (lambda1 ^ 4) * cos(lambda1) * exp((3 * lambda1)) + 0.48e2 * EIy * (lambda1 ^ 4) * sin(lambda1) * exp((3 * lambda1)) - 0.12e2 * EIy * (lambda1 ^ 4) - 0.12e2 * EIy * (lambda1 ^ 4) * (c1 ^ 2) - 0.24e2 * EIy * (lambda1 ^ 4) * c1 + 0.96e2 * EIy * (lambda1 ^ 5) * exp((2 * lambda1)) + 0.12e2 * EIy * (lambda1 ^ 4) * exp((4 * lambda1)) - 0.48e2 * EIy * (lambda1 ^ 4) * (c1 ^ 2) * exp((2 * lambda1)) * cos(lambda1) * sin(lambda1) + 0.32e2 * m * (Omega ^ 2) * (R ^ 4) * (lambda1 ^ 3) * (c1 ^ 2) * exp((2 * lambda1)) - 0.48e2 * EIy * (lambda1 ^ 4) * (c1 ^ 2) * exp((3 * lambda1)) * cos(lambda1) + 0.48e2 * EIy * (lambda1 ^ 4) * exp(lambda1) * (c1 ^ 2) * cos(lambda1) - 0.96e2 * EIy * (lambda1 ^ 4) * c1 * exp((3 * lambda1)) * sin(lambda1) + 0.48e2 * EIy * (lambda1 ^ 4) * (c1 ^ 2) * exp((3 * lambda1)) * sin(lambda1) + 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * c1 * sin(lambda1) * exp((3 * lambda1)) * lambda1 - 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * cos(lambda1) * lambda1 * exp(lambda1) - 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * (lambda1 ^ 2) * c1 * exp((2 * lambda1)) + 0.6e1 * m * (Omega ^ 2) * (R ^ 4) * exp((4 * lambda1)) * (c1 ^ 2) * lambda1 + 0.24e2 * m * (Omega ^ 2) * (R ^ 4) * cos(lambda1) ^ 2 * lambda1 * exp((2 * lambda1)) - 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * cos(lambda1) * exp((3 * lambda1)) * lambda1 - 0.12e2 * m * (Omega ^ 2) * (R ^ 4) * exp((4 * lambda1)) * c1 * lambda1 + 0.12e2 * m * (Omega ^ 2) * (R ^ 4) * (c1 ^ 2) * exp((2 * lambda1)) * lambda1 - 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * c1 * cos(lambda1) * sin(lambda1) * lambda1 * exp((2 * lambda1)) + 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * c1 * sin(lambda1) * lambda1 * exp(lambda1) + 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * (c1 ^ 2) * sin(lambda1) * lambda1 * exp(lambda1) - 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * c1 * cos(lambda1) * lambda1 * exp(lambda1) + 0.6e1 * m * (Omega ^ 2) * (R ^ 4) * lambda1 * exp((4 * lambda1)) + 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * exp((3 * lambda1)) * c1 * cos(lambda1) * lambda1 - 0.24e2 * m * (Omega ^ 2) * (R ^ 4) * (c1 ^ 2) * exp((2 * lambda1)) * lambda1 * cos(lambda1) ^ 2 - 0.48e2 * m * (Omega ^ 2) * (R ^ 4) * (c1 ^ 2) * sin(lambda1) * exp((3 * lambda1)) * lambda1 + 0.60e2 * m * (Omega ^ 2) * (R ^ 4) * exp((2 * lambda1)) * c1 - 0.15e2 * m * (Omega ^ 2) * (R ^ 4) * (c1 ^ 2) * exp((4 * lambda1)) + 0.72e2 * m * (Omega ^ 2) * (R ^ 4) * exp(lambda1) * (c1 ^ 2) * cos(lambda1) + 0.72e2 * m * (Omega ^ 2) * (R ^ 4) * cos(lambda1) * exp((3 * lambda1)) + 0.72e2 * m * (Omega ^ 2) * (R ^ 4) * (c1 ^ 2) * sin(lambda1) * exp(lambda1) + 0.72e2 * m * (Omega ^ 2) * (R ^ 4) * sin(lambda1) * exp((3 * lambda1)) + 0.60e2 * m * (Omega ^ 2) * (R ^ 4) * (c1 ^ 2) * exp((2 * lambda1)) * cos(lambda1) * sin(lambda1) + 0.30e2 * m * (Omega ^ 2) * (R ^ 4) * c1 * exp((4 * lambda1)) - 0.72e2 * m * (Omega ^ 2) * (R ^ 4) * (c1 ^ 2) * exp((3 * lambda1)) * cos(lambda1) - 0.108e3 * m * (Omega ^ 2) * (R ^ 4) * lambda1 * exp((2 * lambda1)) + 0.72e2 * m * (Omega ^ 2) * (R ^ 4) * (c1 ^ 2) * exp((3 * lambda1)) * sin(lambda1) - 0.60e2 * m * (Omega ^ 2) * (R ^ 4) * cos(lambda1) * sin(lambda1) * exp((2 * lambda1)) + 0.144e3 * m * (Omega ^ 2) * (R ^ 4) * c1 * sin(lambda1) * exp(lambda1) - 0.120e3 * m * (Omega ^ 2) * (R ^ 4) * c1 * cos(lambda1) ^ 2 * exp((2 * lambda1)) - 0.72e2 * m * (Omega ^ 2) * (R ^ 4) * cos(lambda1) * exp(lambda1) + 0.72e2 * m * (Omega ^ 2) * (R ^ 4) * exp(lambda1) * sin(lambda1) - 0.144e3 * m * (Omega ^ 2) * (R ^ 4) * c1 * exp((3 * lambda1)) * sin(lambda1) + 0.48e2 * EIy * (lambda1 ^ 4) * cos(lambda1) * sin(lambda1) * exp((2 * lambda1)) + 0.96e2 * EIy * (lambda1 ^ 4) * c1 * cos(lambda1) ^ 2 * exp((2 * lambda1)) + 0.96e2 * EIy * (lambda1 ^ 4) * c1 * sin(lambda1) * exp(lambda1) + 0.48e2 * EIy * (lambda1 ^ 4) * (c1 ^ 2) * sin(lambda1) * exp(lambda1)) / (R ^ 3) / lambda1 / 0.96e2 0
  0 0 (m * Omega ^ 2 * a_cg ^ 2 * R) / 0.2e1 + GK * pi ^ 2 / R / 0.8e1
  ];
% Blade gyroscopic matrix
Gb = [
    0 0 0
    0 0 0
    0 0 0
    ];
% Nacelle/tower top mass matrix
Mt = [
    3*m*R+M, 0, 0, 3*m*R*Ls, 0, 0 
    0, 3*m*R+M, 0, 0, 0, 0
    0, 0, 0.5*m*R*(R^2+3*a_cg^2)+3*m*R*Ls^2+Ix, 0, 0, 0
    3*m*R*Ls, 0, 0, 0.5*m*R*(R^2+3*a_cg^2)+3*m*R*Ls^2+Iz, 0, 0
    0, 0, 0, 0, m*R^3+3*m*R*a_cg^2+Iy, m*R*(R^2+3*a_cg^2)
    0, 0, 0, 0, m*R*(R^2+3*a_cg^2), m*R*(R^2+3*a_cg^2)
    ];
% Nacelle/tower top gyroscopic matrix
Gt = [
    0, 0, 0, 0, 0, 0
    0, 0, 0, 0, 0, 0
    0, 0, 0, -m*R*Omega*(R^2+3*a_cg^2), 0, 0
    0, 0, m*R*Omega*(R^2+3*a_cg^2), 0, 0, 0
    0, 0, 0, 0, 0, 0
    0, 0, 0, 0, 0, 0
    ];
% Nacelle/tower top stiffness matrix
Kt = [
    Kx, 0, 0, 0, gxy, 0
    0, Ky, -gxy, 0, 0, 0
    0, -gxy, Gx, 0, 0, 0
    0, 0, 0, Gz, 0, 0
    gxy, 0, 0, 0, Gy, 0
    0, 0, 0, 0, 0, Gs
    ];