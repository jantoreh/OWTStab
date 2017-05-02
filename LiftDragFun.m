function [Cl,Cd,Clp,Cdp] = LiftDragFun(alpha)
% Aerodynamic paramters
CDfric = 0.005;
alpha_stall = 14*pi/180;
Dalpha_stall = 3*pi/180;
% Normal force coefficient and its derivative
Cn = 2.25*2*pi*sin(alpha)./(4+pi*abs(sin(alpha)));
Cnp = 4.50*pi*cos(alpha)./(4+pi*abs(sin(alpha)))-4.50*pi^2*sin(alpha)./(4+pi*abs(sin(alpha))).^2.*sign(sin(alpha)).*cos(alpha);
% Separation function and its derivative
f = 0.5+0.5*tanh((alpha_stall-abs(alpha))/Dalpha_stall);
fp = -1/2*(1-tanh((alpha_stall-abs(alpha))./Dalpha_stall).^2).*sign(alpha)/Dalpha_stall;
% Aerodynamic coefficients
Cl =  2*pi*alpha.*f + Cn.*cos(alpha).*(1-f);
Cd =  CDfric*f + 3/4*Cn.*sin(alpha).*(1-f);
% Derivatives
Clp = 2*pi*f+2*pi*alpha.*fp+Cnp.*cos(alpha).*(1-f)-Cn.*sin(alpha).*(1-f)-Cn.*cos(alpha).*fp;
Cdp = CDfric*fp+3/4*Cnp.*sin(alpha).*(1-f)+3/4*Cn.*cos(alpha).*(1-f)-3/4*Cn.*sin(alpha).*fp;
