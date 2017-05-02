function AeroMat = AerodynMatrix(par,Omega,W)
% computation parameters
N = 500;
z = linspace(0,par.R,N);
h = z(2)-z(1);

% bending modes
phi_b = BendingModeFun(z,1);
phi_t = TorsionalModeFun(z,1);

Cab=0;Cat=0;Kab=0;Kat=0;Ka1_bt=0;Kb1_bt=0;Ka0_tb=0;Ka1_tb=0;Kb1_tb=0;
Ca0_tb=0;Ca1_tb=0;Cb1_tb=0;Ca0_bt=0;Ca1_bt=0;Cb1_bt=0;
for i=1:N
    
    dAeroMat = dAerodynMatrix(par,z(i),Omega,W,phi_b(i),phi_t(i));
    if i > 2
        Cab      = Cab    + h/2*(dAeroMat_old.Cab    + dAeroMat.Cab    );
        Cat      = Cat    + h/2*(dAeroMat_old.Cat    + dAeroMat.Cat    );
        Kab      = Kab    + h/2*(dAeroMat_old.Kab    + dAeroMat.Kab    );
        Kat      = Kat    + h/2*(dAeroMat_old.Kat    + dAeroMat.Kat    );
        Ka1_bt    = Ka1_bt  + h/2*(dAeroMat_old.Ka1_bt  + dAeroMat.Ka1_bt  );
        Kb1_bt    = Kb1_bt  + h/2*(dAeroMat_old.Kb1_bt  + dAeroMat.Kb1_bt  );
        Ka0_tb   = Ka0_tb + h/2*(dAeroMat_old.Ka0_tb + dAeroMat.Ka0_tb );
        Ka1_tb   = Ka1_tb + h/2*(dAeroMat_old.Ka1_tb + dAeroMat.Ka1_tb );
        Kb1_tb   = Kb1_tb + h/2*(dAeroMat_old.Kb1_tb + dAeroMat.Kb1_tb );
        Ca0_tb   = Ca0_tb + h/2*(dAeroMat_old.Ca0_tb + dAeroMat.Ca0_tb );
        Ca1_tb   = Ca1_tb + h/2*(dAeroMat_old.Ca1_tb + dAeroMat.Ca1_tb );
        Cb1_tb   = Cb1_tb + h/2*(dAeroMat_old.Cb1_tb + dAeroMat.Cb1_tb );
        Ca0_bt   = Ca0_bt + h/2*(dAeroMat_old.Ca0_bt + dAeroMat.Ca0_bt );
        Ca1_bt   = Ca1_bt + h/2*(dAeroMat_old.Ca1_bt + dAeroMat.Ca1_bt );
        Cb1_bt   = Cb1_bt + h/2*(dAeroMat_old.Cb1_bt + dAeroMat.Cb1_bt );
    end
    dAeroMat_old = dAeroMat;
end

AeroMat.Cab    = Cab   ;
AeroMat.Kab    = Kab   ;
AeroMat.Cat    = Cat   ;
AeroMat.Kat    = Kat   ;
AeroMat.Ka1_bt = Ka1_bt ;
AeroMat.Kb1_bt = Kb1_bt ;
AeroMat.Ka0_tb = Ka0_tb;
AeroMat.Ka1_tb = Ka1_tb;
AeroMat.Kb1_tb = Kb1_tb;
AeroMat.Ca0_tb = Ca0_tb;
AeroMat.Ca1_tb = Ca1_tb;
AeroMat.Cb1_tb = Cb1_tb;
AeroMat.Ca0_bt = Ca0_bt;
AeroMat.Ca1_bt = Ca1_bt;
AeroMat.Cb1_bt = Cb1_bt;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dAeroMat = dAerodynMatrix(par,z,Omega,W,phi_b,phi_t)
% aerodynamic paramters
c = par.c;
Ls = par.Ls;
rho = par.rho;

% computer parameters
lambda = z*Omega/W;
alpha0 = atan(1/lambda);
psi0 = alpha0; % zero twist of blade
U0 = sqrt(W^2+(z*Omega)^2);

[Cl,Cd,dCl,dCd] = LiftDragFun(alpha0);



CX0 = Cl*sin(psi0)-Cd*cos(psi0);
CY0 = Cl*cos(psi0)+Cd*sin(psi0);
CXp0 = dCl*sin(psi0)-dCd*cos(psi0);
CYp0 = dCl*cos(psi0)+dCd*sin(psi0);


% computer matrix
dAeroMat.Cab = 0.25*c*rho*W*[
    4*phi_b^2*CY0+2*phi_b^2*CYp0*lambda-2*phi_b^2*CX0*lambda, 2*phi_b*phi_b*CYp0-2*phi_b*phi_b*CX0-4*phi_b*phi_b*CY0*lambda, -c*phi_b*phi_t*CYp0*lambda
    4*phi_b*phi_b*CX0+2*phi_b*phi_b*CXp0*lambda+2*phi_b*phi_b*CY0*lambda, 2*phi_b^2*CXp0+2*phi_b^2*CY0-4*phi_b^2*CX0*lambda, -c*phi_b*phi_t*CXp0*lambda
    0, 0, 0
    ];

dAeroMat.Kab = 0.5*c*rho*U0^2*[
    0, 0, -phi_t*phi_b*CYp0
    0, 0, -phi_b*phi_t*CXp0
    0, 0, 0
    ];



dAeroMat.Cat = 0.125*c*rho*W*[
    6*CXp0+6*CY0-12*CX0*lambda, 0, -6*lambda*z*CY0-6*lambda*z*CXp0-12*z*CX0, 6*Ls*CY0+(6*Ls-3*c*lambda)*CXp0-12*Ls*CX0*lambda, 0, 0
    0, 24*CY0+12*CYp0*lambda-12*CX0*lambda, 0, 0, -24*lambda*z*CY0+12*z*CYp0-12*z*CX0, -24*lambda*z*CY0+12*z*CYp0-12*z*CX0
    12*lambda*z*CY0-6*z*CYp0+6*z*CX0, 0, (6*Ls^2+12*z^2)*CY0+6*z^2*CYp0*lambda+(6*Ls^2-3*Ls*c*lambda)*CXp0+(-6*z^2*lambda-12*Ls^2*lambda)*CX0, 18*Ls*lambda*z*CY0+3*(-2*Ls+c*lambda)*z*CYp0+6*Ls*lambda*z*CXp0+18*z*Ls*CX0, 0, 0
    6*Ls*CXp0+6*Ls*CY0-12*Ls*CX0*lambda, 0, -18*Ls*lambda*z*CY0-3*(-2*Ls+c*lambda)*z*CYp0-6*Ls*lambda*z*CXp0-18*z*Ls*CX0, (6*Ls^2+12*z^2)*CY0+6*z^2*CYp0*lambda+(6*Ls^2-3*Ls*c*lambda)*CXp0+(-6*z^2*lambda-12*Ls^2*lambda)*CX0, 0, 0
    0, 12*lambda*z*CY0+12*lambda*z*CXp0+24*z*CX0, 0, 0, 12*z^2*CY0+12*z^2*CXp0-24*z^2*CX0*lambda, 12*z^2*CY0+12*z^2*CXp0-24*z^2*CX0*lambda
    0, 12*lambda*z*CY0+12*lambda*z*CXp0+24*z*CX0, 0, 0, 12*z^2*CY0+12*z^2*CXp0-24*z^2*CX0*lambda, 12*z^2*CY0+12*z^2*CXp0-24*z^2*CX0*lambda
    ];



dAeroMat.Kat = 0.25*c*rho*W*W*[
    0, 0, 0, 3*CY0+6*CY0*lambda^2+6*CX0*lambda-3*CXp0, 0, 0
    0, 0, 0, 0, 0, 0
    0, 0, -3*Ls*CY0+6*Ls*CX0*lambda-3*Ls*CXp0, (6*CX0*lambda^2-6*CY0*lambda+3*CYp0+3*CX0)*z, 0, 0
    0, 0, (-3*CYp0+6*CY0*lambda+3*CX0)*z, -3*Ls*CY0+6*Ls*CX0*lambda-3*Ls*CXp0, 0, 0
    0, 0, 0, 0, 0, 0
    0, 0, 0, 0, 0, 0
    ];

dAeroMat.Ka1_bt = 0.5*c*rho*W*W*[
    0, 0, 0, -phi_b*CYp0+2*phi_b*CY0*lambda+phi_b*CX0, 0, 0
    0, 0, 0, -phi_b*CY0+2*phi_b*CX0*lambda-phi_b*CXp0, 0, 0
    0, 0, 0, 0, 0, 0
    ];

dAeroMat.Kb1_bt = 0.5*c*rho*W*W*[
    0, 0, -phi_b*CYp0+2*phi_b*CY0*lambda+phi_b*CX0, 0, 0, 0
    0, 0, -phi_b*CY0+2*phi_b*CX0*lambda-phi_b*CXp0, 0, 0, 0
    0, 0, 0, 0, 0, 0
    ];

dAeroMat.Ka0_tb = 0.5*c*rho*U0^2*[
    0, 0, 0
    0, 0, -phi_t*CYp0
    0, 0, 0
    0, 0, 0
    0, 0, -z*phi_t*CXp0
    0, 0, -z*phi_t*CXp0
    ];

dAeroMat.Ka1_tb = 0.5*c*rho*U0^2*[
    0, 0, -phi_t*CXp0
    0, 0, 0
    0, 0, phi_t*z*CYp0
    phi_b*CX0, -phi_b*CY0, -phi_t*Ls*CXp0
    0, 0, 0
    0, 0, 0
    ];

dAeroMat.Kb1_tb = 0.5*c*rho*U0^2*[
    0, 0, 0
    0, 0, 0
    phi_b*CX0, -phi_b*CY0, -phi_t*Ls*CXp0
    0, 0, -phi_t*z*CYp0
    0, 0, 0
    0, 0, 0
    ];

dAeroMat.Ca0_tb = 0.5*c*rho*W*[
    0, 0, 0
    2*phi_b*CY0+phi_b*CYp0*lambda-phi_b*CX0*lambda, phi_b*CYp0-phi_b*CX0-2*phi_b*CY0*lambda, -1/2*c*phi_t*CYp0*lambda
    0, 0, 0
    0, 0, 0
    phi_b*(2*CX0+CXp0*lambda+CY0*lambda)*z, -phi_b*(-CXp0-CY0+2*CX0*lambda)*z, -1/2*c*phi_t*CXp0*lambda*z
    phi_b*(2*CX0+CXp0*lambda+CY0*lambda)*z, -phi_b*(-CXp0-CY0+2*CX0*lambda)*z, -1/2*c*phi_t*CXp0*lambda*z
    ];

dAeroMat.Ca1_tb = 0.5*c*rho*W*[
    2*phi_b*CX0+phi_b*CXp0*lambda+phi_b*CY0*lambda, phi_b*CXp0+phi_b*CY0-2*phi_b*CX0*lambda, -1/2*c*phi_t*CXp0*lambda
    0, 0, 0
    -phi_b*(2*CY0+CYp0*lambda-CX0*lambda)*z, phi_b*(-CYp0+CX0+2*CY0*lambda)*z, 1/2*c*phi_t*CYp0*lambda*z
    2*phi_b*Ls*CX0+phi_b*Ls*CXp0*lambda+phi_b*Ls*CY0*lambda, phi_b*Ls*CXp0+phi_b*Ls*CY0-2*phi_b*Ls*CX0*lambda, -1/2*c*phi_t*Ls*CXp0*lambda
    0, 0, 0
    0, 0, 0
    ];

dAeroMat.Cb1_tb = 0.5*c*rho*W*[
    0, 0, 0
    0, 0, 0
    2*phi_b*Ls*CX0+phi_b*Ls*CXp0*lambda+phi_b*Ls*CY0*lambda, phi_b*Ls*CXp0+phi_b*Ls*CY0-2*phi_b*Ls*CX0*lambda, -1/2*c*phi_t*Ls*CXp0*lambda
    phi_b*(2*CY0+CYp0*lambda-CX0*lambda)*z, -phi_b*(-CYp0+CX0+2*CY0*lambda)*z, -1/2*c*phi_t*CYp0*lambda*z
    0, 0, 0
    0, 0, 0
    ];

dAeroMat.Ca0_bt = 0.5*c*rho*W*[
    0, 2*phi_b*CY0+phi_b*CYp0*lambda-phi_b*CX0*lambda, 0, 0, -phi_b*(-CYp0+CX0+2*CY0*lambda)*z, -phi_b*(-CYp0+CX0+2*CY0*lambda)*z
    0, 2*phi_b*CX0+phi_b*CXp0*lambda+phi_b*CY0*lambda, 0, 0, -phi_b*(-CXp0-CY0+2*CX0*lambda)*z, -phi_b*(-CXp0-CY0+2*CX0*lambda)*z
    0, 0, 0, 0, 0, 0
    ];

dAeroMat.Ca1_bt = 0.5*c*rho*W*[
    phi_b*CYp0-phi_b*CX0-2*phi_b*CY0*lambda, 0, -phi_b*(2*CY0+CYp0*lambda-CX0*lambda)*z, -phi_b*Ls*CX0+phi_b*Ls*CYp0-1/2*phi_b*c*CYp0*lambda-2*phi_b*Ls*CY0*lambda, 0, 0
    phi_b*CXp0+phi_b*CY0-2*phi_b*CX0*lambda, 0, -phi_b*(2*CX0+CXp0*lambda+CY0*lambda)*z, phi_b*Ls*CXp0+phi_b*Ls*CY0-1/2*phi_b*c*CXp0*lambda-2*phi_b*Ls*CX0*lambda, 0, 0
    0, 0, 0, 0, 0, 0
    ];

dAeroMat.Cb1_bt = 0.5*c*rho*W*[
    0, 0, -phi_b*Ls*CX0+phi_b*Ls*CYp0-1/2*phi_b*c*CYp0*lambda-2*phi_b*Ls*CY0*lambda, phi_b*(2*CY0+CYp0*lambda-CX0*lambda)*z, 0, 0
    0, 0, phi_b*Ls*CXp0+phi_b*Ls*CY0-1/2*phi_b*c*CXp0*lambda-2*phi_b*Ls*CX0*lambda, phi_b*(2*CX0+CXp0*lambda+CY0*lambda)*z, 0, 0
    0, 0, 0, 0, 0, 0
    ];







