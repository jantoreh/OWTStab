
% DTU course main file

% Team magenta: Jan-Tore & Joey

clc
clear
close all

% ------------------------------------------------ %
% Define parameters
par        = Parameters();
Omega      = 1.5;
Wsp        = 5:25;
nfreq      = 11;
AeroSwitch = 1; % Aerodynamics on/off

% ------------------------------------------------ %
% Eigen value analysis as function of omega

fid=fopen('modes.dat','w');
for k = 1:length(Omega)
for i = 1:length(Wsp)

    % Aerodynamics
    AeroMat = AerodynMatrix(par,Omega(k),Wsp(i));
    
    % Find matrices for blades
    [Mb,Kb,Gb,Mt,Gt,Kt]=DiagonalStructuralMatrices(par,Omega(k));
    ndofb = size(Mb,1);
    
    % Find matrices for tower top
    psi0 = 0;
    psi1 = psi0;
    psi2 = psi0 + 2*pi/3;
    psi3 = psi0 + 2*2*pi/3;
    [Mtb1,Gbt1,Gtb1,Ktb1]=AzimutDependentStructuralMatrix(par,Omega(k),psi1);
    [Mtb2,Gbt2,Gtb2,Ktb2]=AzimutDependentStructuralMatrix(par,Omega(k),psi2);
    [Mtb3,Gbt3,Gtb3,Ktb3]=AzimutDependentStructuralMatrix(par,Omega(k),psi3);
    ndoft = size(Mtb1,1);
    ndof = 3*ndofb + ndoft; % Three blades + tower top

    % Assemble total matrices
    z = zeros(ndofb);
    z2=zeros(ndofb,ndoft);
    I = eye(3);
    
    % Azimuth dependent aerodynamic matrices
    Kab = AeroMat.Kab*AeroSwitch;
    Cab = AeroMat.Cab*AeroSwitch;
    Cat = AeroMat.Cat*AeroSwitch;
    Kat = AeroMat.Kat*AeroSwitch;
    Ca0_bt=AeroMat.Ca0_bt*AeroSwitch;
    Ca1_bt=AeroMat.Ca1_bt*AeroSwitch;
    Cb1_bt=AeroMat.Cb1_bt*AeroSwitch;
    Ca0_tb=AeroMat.Ca0_tb*AeroSwitch;
    Ca1_tb=AeroMat.Ca1_tb*AeroSwitch;
    Cb1_tb=AeroMat.Cb1_tb*AeroSwitch;
    Ka1_bt=AeroMat.Ka1_bt*AeroSwitch;
    Kb1_bt=AeroMat.Kb1_bt*AeroSwitch;
    Ka0_tb=AeroMat.Ka0_tb*AeroSwitch;
    Ka1_tb=AeroMat.Ka1_tb*AeroSwitch;
    Kb1_tb=AeroMat.Kb1_tb*AeroSwitch;
    
    Kabt1 = Ka1_bt*cos(psi1) + Kb1_bt*sin(psi1);
    Kabt2 = Ka1_bt*cos(psi2) + Kb1_bt*sin(psi2);
    Kabt3 = Ka1_bt*cos(psi3) + Kb1_bt*sin(psi3);
    
    Katb1 = Ka0_tb + Ka1_tb*cos(psi1) + Kb1_tb*sin(psi1);
    Katb2 = Ka0_tb + Ka1_tb*cos(psi2) + Kb1_tb*sin(psi2);
    Katb3 = Ka0_tb + Ka1_tb*cos(psi3) + Kb1_tb*sin(psi3);
    
    Cabt1 = Ca0_bt + Ca1_bt*cos(psi1) + Cb1_bt*sin(psi1);
    Cabt2 = Ca0_bt + Ca1_bt*cos(psi2) + Cb1_bt*sin(psi2);
    Cabt3 = Ca0_bt + Ca1_bt*cos(psi3) + Cb1_bt*sin(psi3);
    
    Catb1 = Ca0_tb + Ca1_tb*cos(psi1) + Cb1_tb*sin(psi1);
    Catb2 = Ca0_tb + Ca1_tb*cos(psi2) + Cb1_tb*sin(psi2);
    Catb3 = Ca0_tb + Ca1_tb*cos(psi3) + Cb1_tb*sin(psi3);
    
    % Construct system matrices
    M = [ Mb,   z,   z, Mtb1';...
           z,  Mb,   z, Mtb2';...
           z,   z,  Mb, Mtb3';...
        Mtb1,Mtb2,Mtb3,   Mt];
    
    C = [Gb+Cab,         z,          z, Gbt1+Cabt1;...
              z,    Gb+Cab,          z, Gbt2+Cabt2;...
              z,         z,     Gb+Cab, Gbt3+Cabt3;...
     Gtb1+Catb1,Gtb2+Catb2, Gtb3+Catb3,    Gt+Cat];
       
    K = [Kb+Kab,          z,          z, Kabt1;...
              z ,    Kb+Kab,          z, Kabt2;...
              z,          z,     Kb+Kab, Kabt3;...
     Ktb1+Katb1, Ktb2+Katb2, Ktb3+Katb3, Kt+Kat];
    
    
    
    B = [ I, I.*cos(psi1), I.*sin(psi1),     z2;...
          I, I.*cos(psi2), I.*sin(psi2),     z2;...
          I, I.*cos(psi3), I.*sin(psi3),     z2;...
        z2',          z2',          z2',eye(6)];
    
    mu = [1/3*I,    z,     z,  z2;...
              z,2/3*I,     z,  z2;...
              z,    z, 2/3*I,  z2;...
            z2',  z2',   z2', eye(6)];
    
    R = [zeros(3,15);...
         z,z,Omega*I,z2;...
         z,-Omega*I,z,z2;...
         zeros(6,15)];
    
    KB = mu*B'*K*B;
    CB = mu*B'*C*B;
    MB = mu*B'*M*B;
    
    
    % Construct AB
    AB = [zeros(ndof),eye(ndof);...
          -MB\(MB*R.^2+CB*R+KB),-MB\(2*MB*R+CB)];
    
    % Eigenvalue analysis
    [Vb,Db]=eig(AB);
    
    
    
    % Eigen values
    lambda = diag(Db);
    [tmp,I] = sort(imag(lambda));
    lambda=lambda(I);

    id = tmp>0;
    
    lambda    = lambda(id);
    tmp       = tmp(id)/2/pi; % Pick positive frequencies
    Freq(i,:,k) = tmp(1:nfreq);
    tmp       = -real(lambda)./abs(lambda);
    Damp(i,:,k) = tmp(1:nfreq);
    
    
    % Eigen vectors
    Vb             = Vb(:,I); % Sort vectors according to ascending frequencies
    Vec            = Vb(1:ndof,(ndof+1):(ndof+nfreq)); % Pick vectors for wanted frequencies
    ModeVec(:,:,i,k) = Vec; % Store for later use
    
    
    fprintf(fid,' %f ',Wsp(i));
    for ii=1:size(Vec,2)
        w0         = Vec(1:3,ii);
        wc         = Vec(4:6,ii);
        ws         = Vec(7:9,ii);
        A0_amp     = abs(w0);
        ABW_amp    = 0.5*abs(wc-sqrt(-1)*ws);
        AFW_amp    = 0.5*abs(wc+sqrt(-1)*ws);
        A0_arg     = angle(w0);
        ABW_arg    = angle(wc-sqrt(-1)*ws);
        AFW_arg    = angle(wc+sqrt(-1)*ws);
        g_amp      = abs(Vec(10:15,ii));
        g_arg      = angle(Vec(10:15,ii));
        
        % Rotor components
        for j=1:3
          fprintf(fid,' %e %e %e %e %e %e ',A0_amp(j),A0_arg(j),ABW_amp(j),ABW_arg(j),AFW_amp(j),AFW_arg(j));
        end
        % Ground 
        for j=1:length(g_amp)
            fprintf(fid,' %e %e ',g_amp(j),g_arg(j));
        end

    end
    fprintf(fid,'\n');
    
    
    
end
end
fclose(fid);



% Save frequencies and damping - ascending
tmp = [];
for i=1:nfreq
    tmp = [tmp,Freq(:,i,1),Damp(:,i,1)];
end
tmp = [Wsp',tmp];
save('freq.dat','tmp','-ascii')

% Load structural frequencies
f = load('freq_struct.dat');
f = f(end,2:end);
f = repmat(f,size(Freq,1),1);

% Plot
labels={'Tower side-side','Tower fore-aft','Sym. edge/DT','BW flapwise',...
    'Sym. flap','FW flapwise','BW edgewise','FW edgewise','BW flapwise 2',...
    'Sym. flap 2','FW flapwise 2'};
figure
subplot(1,2,1)
color=distinguishable_colors(11);
for i=1:nfreq
    plot(Wsp,Freq(:,i),'-o','color',color(i,:));hold on
end
for i=1:nfreq
    plot(Wsp,f(:,i),'-x','color',color(i,:));hold on
end
legend(labels)
xlabel('Wind speed [m/s]')
ylabel('Frequency [Hz]')

subplot(1,2,2)
for i=1:nfreq
    plot(Wsp,Damp(:,i),'-o','color',color(i,:));hold on
end
legend(labels)
xlabel('Wind speed [rad/s]')
ylabel('Damping ratio [-]')

stop

% ------------------------------------------------ %
%% Save results for delivery

m = load('modes.dat');

% Normalize
%for i=1:nfreq
%    id = 1:2:29;
%    m(:,30*i+1-30+id) = m(:,30*i+1-30+id)./max(max(m(:,30*i+1-30+id)));
%end
%save('modes_norm.dat','m','-ascii');


modeno = 3;
col1   = modeno*30-28;

% Plot wanted mode
o = m(:,1); % Omega
figure('name',['Mode ' num2str(modeno)])
subplot(3,2,1);hold on
id = 1;
plot(o,m(:,col1-1+id),'-ok')
plot(o,m(:,col1-1+id+2),'-xk')
plot(o,m(:,col1-1+id+4),'-dk')
legend('Sym.','Back.','Forw.')
title('Flapwise')
ylabel('Amplitude')
yyaxis right
id = id+1;
plot(o,m(:,col1-1+id),'-o')
plot(o,m(:,col1-1+id+2),'-x')
plot(o,m(:,col1-1+id+4),'-d')
ylabel('Phase')

subplot(3,2,2);hold on
id = 7;
plot(o,m(:,col1-1+id),'-ok')
plot(o,m(:,col1-1+id+2),'-xk')
plot(o,m(:,col1-1+id+4),'-dk')
legend('Sym.','Back.','Forw.')
title('Edgewise')
ylabel('Amplitude')
yyaxis right
id = id+1;
plot(o,m(:,col1-1+id),'-o')
plot(o,m(:,col1-1+id+2),'-x')
plot(o,m(:,col1-1+id+4),'-d')
ylabel('Phase')

subplot(3,2,3);hold on
id = 13;
plot(o,m(:,col1-1+id),'-ok')
plot(o,m(:,col1-1+id+2),'-xk')
plot(o,m(:,col1-1+id+4),'-dk')
legend('Sym.','Back.','Forw.')
title('Torsion')
ylabel('Amplitude')
yyaxis right
id = id+1;
plot(o,m(:,col1-1+id),'-o')
plot(o,m(:,col1-1+id+2),'-x')
plot(o,m(:,col1-1+id+4),'-d')
ylabel('Phase')

subplot(3,2,4);hold on
id = 19;
plot(o,m(:,col1-1+id),'-ok')
plot(o,m(:,col1-1+id+2),'-xk')
title('Tower/nacelle trans.')
legend('Lat.','Long.')
ylabel('Amplitude')
yyaxis right
id = id+1;
plot(o,m(:,col1-1+id),'-o')
plot(o,m(:,col1-1+id+2),'-x')
ylabel('Phase')

subplot(3,2,5);hold on
id = 23;
plot(o,m(:,col1-1+id),'-ok')
plot(o,m(:,col1-1+id+2),'-xk')
title('Nacelle tilt/yaw')
legend('Tilt','Yaw')
ylabel('Amplitude')
yyaxis right
id = id+1;
plot(o,m(:,col1-1+id),'-o')
plot(o,m(:,col1-1+id+2),'-x')
ylabel('Phase')

subplot(3,2,6);hold on
id = 27;
plot(o,m(:,col1-1+id),'-ok')
plot(o,m(:,col1-1+id+2),'-xk')
title('Nacelle roll/ DT torsion')
legend('Roll','Tors.')
ylabel('Amplitude')
yyaxis right
id = id+1;
plot(o,m(:,col1-1+id),'-o')
plot(o,m(:,col1-1+id+2),'-x')
ylabel('Phase')



% ------------------------------------------------ %
%% Animation

modeno = 1;
omegano = 1;
wspno = 10;


% Mode shapes
N = 50; % Number of points per mode shape
r = linspace(0,50,N); % Blades
% Blade bending
phib = BendingModeFun(r,1);
% Blade torsion
phit = TorsionalModeFun(r,1);


% Amplitudes
ModeAmp = abs(ModeVec(:,modeno,wspno));
% Phases
ModePhase = angle(ModeVec(:,modeno,wspno));

% Amplitude scale
AmpScale = 50;

[xb,yb,zb] = BladeGen(r);
xb=-xb+1;
yb=yb*2;

%xb = zeros(size(r));
%yb = zeros(size(r));
%zb = r;

TT = 70; % Tower top
Overhang = 1;


t = 0;
dt = 0.1;
c = 3; % Chord length
xb = xb*2;

figure('units','normalized','outerposition',[0 0 1 1]);
ii = 1;
while t < 1000
    t = t+dt;
    
    for i=1:length(modeno)
        
    % Amplitudes
    ModeAmp = abs(ModeVec(:,modeno(i),wspno));
    % Phases
    ModePhase = angle(ModeVec(:,modeno(i),wspno));
    
    % Nacelle motion
    ss = -AmpScale*ModeAmp(10,:).*0.5*phib.*cos(Omega(omegano)*t+ModePhase(10,:));
    fa = -AmpScale*ModeAmp(11,:).*0.5*phib.*cos(Omega(omegano)*t+ModePhase(11,:));
    % Global coordinates
    zt = r*TT/max(r);
    xt = ss;
    yt = fa;
    
    % Blade 1
    f1 = AmpScale*ModeAmp(1,:).*phib.*cos(Omega(omegano)*t+ModePhase(1,:));
    e1 = AmpScale*ModeAmp(2,:).*phib.*cos(Omega(omegano)*t+ModePhase(2,:));
    rot1 = AmpScale*ModeAmp(3,:).*phit.*cos(Omega(omegano)*t+ModePhase(3,:));
    % In local coordinates
    xb1l = xb + e1;
    yb1l = yb + f1 + rot1.*xb.*c.*(xb-0.5);
    zb1l = zb;
    % In global coordinates
    xb1g = xb1l.*cos(psi1) + zb1l.*sin(psi1) + ss(end);
    yb1g = yb1l*cos(psi1) - Overhang + fa(end);
    zb1g = zb1l*cos(psi1) - xb1l*sin(psi1) + TT;
    
    
    % Blade 2
    f2 = AmpScale*ModeAmp(4,:).*phib.*cos(Omega(omegano)*t+ModePhase(4,:));
    e2 = AmpScale*ModeAmp(5,:).*phib.*cos(Omega(omegano)*t+ModePhase(5,:));
    rot2 = AmpScale*ModeAmp(6,:).*phit.*cos(Omega(omegano)*t+ModePhase(6,:));
    % In local coordinates
    xb2l = xb + e2;
    yb2l = yb + f2 + rot2.*xb.*c.*(xb-0.5);
    zb2l = zb;
    % In global coordinates
    xb2g = xb2l.*cos(psi2) + zb2l.*sin(psi2) +ss(end);
    yb2g = yb2l - Overhang + fa(end);
    zb2g = zb2l*cos(psi2) - xb2l*sin(psi2) + TT;
    
    % Blade 2
    f3 = AmpScale*ModeAmp(7,:).*phib.*cos(Omega(omegano)*t+ModePhase(7,:));
    e3 = AmpScale*ModeAmp(8,:).*phib.*cos(Omega(omegano)*t+ModePhase(8,:));
    rot3 = AmpScale*ModeAmp(9,:).*phit.*cos(Omega(omegano)*t+ModePhase(9,:));
    % In local coordinates
    xb3l = xb + e3;
    yb3l = yb + f3 + rot3.*xb.*c.*(xb-0.5);
    zb3l = zb;
    % In global coordinates
    xb3g = xb3l*cos(psi3) + zb3l.*sin(psi3) + ss(end);
    yb3g = yb3l - Overhang + fa(end);
    zb3g = zb3l*cos(psi3) - xb3l*sin(psi3) + TT;
    
    
    
    
    Nr = max([floor(length(modeno)/2),1]);
    Nc = ceil(length(modeno)/Nr);
    subplot(Nr,Nc,i)
    %scatter3(xb1g,yb1g,zb1g,'.k')
    surf(xb1g,yb1g,zb1g,zeros(size(zb)))
    title(['Mode ',num2str(modeno(i))])
    hold on
    %scatter3(xb2g,yb2g,zb2g,'.k')
    %scatter3(xb3g,yb3g,zb3g,'.k')
    surf(xb2g,yb2g,zb2g,zeros(size(zb)))
    surf(xb3g,yb3g,zb3g,zeros(size(zb)))
    scatter3(ss,fa,zt,'.k')
    hold off
    axis([-70,70,-50,50,0,70+TT])
    view([50,30])
    
    
    
    end
    
    pause(0.001)
    
    Movie(ii) = getframe;
    ii=ii+1;
    
    
    
end



