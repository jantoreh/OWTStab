function [Mtb,Gbt,Gtb,Ktb] = AzimutDependentStructuralMatrix(par,Omega,psi)
% Shape functions used in the modal expansion
N = 500;
R = par.R;
z = linspace(0,R,N);
h = z(2)-z(1);
phi_b = BendingModeFun(z,1);
phi_t = TorsionalModeFun(z,1);
% Approximate integrals by trapez formula 
Mtb=0; Gbt=0; Gtb=0; Ktb=0;
[dMtb_old,dGbt_old,dGtb_old,dKtb_old] = dMatrix(z(1),phi_b(1),phi_t(1),psi,par,Omega);
for i=2:N
    [dMtb,dGbt,dGtb,dKtb] = dMatrix(z(i),phi_b(i),phi_t(i),psi,par,Omega);
    Mtb = Mtb + h/2*(dMtb_old+dMtb);
    Gtb = Gtb + h/2*(dGtb_old+dGtb);
    Gbt = Gbt + h/2*(dGbt_old+dGbt);
    Ktb = Ktb + h/2*(dKtb_old+dKtb);
    dMtb_old=dMtb;dGbt_old=dGbt;dGtb_old=dGtb;dKtb_old=dKtb;
end
%==========================================================================
% Kernel of integrals in Mtb, Gbt, Gtb, and Ktb
%==========================================================================
function [Mtb,Gbt,Gtb,Ktb] = dMatrix(z,phi_b,phi_t,psi,par,Omega)
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
% Mass coupling between rotor and tower
Mtb = m*[
    0,cos(psi)*phi_b,0
    phi_b,0,-a_cg*phi_t
    -phi_b*(a_cg*sin(psi)+z*cos(psi)), Ls*sin(psi)*phi_b, a_cg*phi_t*(a_cg*sin(psi)+z*cos(psi))
    -phi_b*(a_cg*cos(psi)-z*sin(psi)), Ls*cos(psi)*phi_b, a_cg*phi_t*(a_cg*cos(psi)-z*sin(psi))
    0, phi_b*z, 0
    0, phi_b*z, 0
    ];
% Gyroscopic coupling from tower to rotor 
Gbt = m*Omega*[
    0, 0, -2*phi_b*(a_cg*cos(psi)-z*sin(psi)), 2*phi_b*(a_cg*sin(psi)+z*cos(psi)), 0, 0
    0, 0, 0, 0, 2*a_cg*phi_b, 2*a_cg*phi_b
    0, 0, 2*a_cg*phi_t*(a_cg*cos(psi)-z*sin(psi)), -2*a_cg*phi_t*(a_cg*sin(psi)+z*cos(psi)), 0, 0
    ];
% Gyroscopic coupling from rotor to tower
Gtb = m*Omega*[
    0, -2*sin(psi)*phi_b, 0
    0, 0, 0
    0, 2*Ls*cos(psi)*phi_b, 0
    0, -2*Ls*sin(psi)*phi_b, 0
    0, -2*a_cg*phi_b, 0
    0, -2*a_cg*phi_b, 0
    ];
% Deflection proportional inertia coupling from rotor to tower
Ktb = m*Omega^2*[
    0, -cos(psi)*phi_b, 0
    0, 0, 0
    -phi_b*(a_cg*sin(psi)+z*cos(psi)), -Ls*sin(psi)*phi_b, a_cg*phi_t*(a_cg*sin(psi)+z*cos(psi))
    -phi_b*(a_cg*cos(psi)-z*sin(psi)), -Ls*cos(psi)*phi_b, a_cg*phi_t*(a_cg*cos(psi)-z*sin(psi))
    0, 0, 0
    0, 0, 0
    ];