
% DTU course main file

% Team magenta: Jan-Tore & Joey

clc
clear
close all

% ------------------------------------------------ %
% Define parameters
par = Parameters();
Omega = 0:0.05:1.5;
r = 0:1:par.R;
nfreq = 11;

% ------------------------------------------------ %
% Eigen value analysis as function of omega

fid=fopen('modes.dat','w');
for i = 1:length(Omega)

    % Find matrices for blades
    [Mb,Kb,Gb,Mt,Gt,Kt]=DiagonalStructuralMatrices(par,Omega(i));
    ndofb = size(Mb,1);
    
    % Find matrices for tower top
    psi0 = 0;
    psi1 = psi0;
    psi2 = psi0 + 2*pi/3;
    psi3 = psi0 + 2*2*pi/3;
    [Mtb1,Gbt1,Gtb1,Ktb1]=AzimutDependentStructuralMatrix(par,Omega(i),psi1);
    [Mtb2,Gbt2,Gtb2,Ktb2]=AzimutDependentStructuralMatrix(par,Omega(i),psi2);
    [Mtb3,Gbt3,Gtb3,Ktb3]=AzimutDependentStructuralMatrix(par,Omega(i),psi3);
    ndof = 3*ndofb + size(Mtb1,1); % Three blades + tower top

    % Assemble total matrices
    z = zeros(ndofb);
    M = [Mb,z,z,Mtb1';z,Mb,z,Mtb2';z,z,Mb,Mtb3';Mtb1,Mtb2,Mtb3,Mt];
    C = [Gb,z,z,Gbt1;z,Gb,z,Gbt2;z,z,z,Gbt3;Gtb1,Gtb2,Gtb3,Gt];
    K = [Kb,z,z,zeros(3,6);z,Kb,z,zeros(3,6);z,z,Kb,zeros(3,6);Ktb1,Ktb2,Ktb3,Kt];
    
    I = eye(3);
    B = [I,I.*cos(psi1),I.*sin(psi1),zeros(3,6);...
        I,I.*cos(psi2),I.*sin(psi2),zeros(3,6);...
        I,I.*cos(psi3),I.*sin(psi3),zeros(3,6);...
        zeros(6,9),eye(6)];
    
    mu = [1/3*I,zeros(3,12);...
          z,2/3*I,zeros(3,9);...
          z,z,2/3*I,zeros(3,6);...
          zeros(6,9),eye(6)];
    
    R = [zeros(3,15);...
         z,z,Omega(i)*I,zeros(3,6);...
         z,-Omega(i)*I,z,zeros(3,6);...
         zeros(6,15)];
    
    KB = mu*B'*K*B;
    CB = mu*B'*C*B;
    MB = mu*B'*M*B;
    
    
    
    AB = [zeros(ndof),eye(ndof);...
          -MB\(MB*R.^2+CB*R+KB),-MB\(2*MB*R+CB)];
    
    % Eigenvalue analysis
    [Vb,Db]=eig(AB);
    
    
    
    % Eigen values
    [tmp,I] = sort(diag(imag(Db))/2/pi); % Sort for ascending frequencies
    Freq(i,:)=tmp(tmp>0); % Pick positive frequencies
    
    
    
    % Eigen vectors
    Vb = Vb(:,I); % Sort vectors according to ascending frequencies
    Vec = Vb(1:ndof,(ndof+1):(2*ndof)); % Pick vectors
    ModeVec(:,:,i) = Vec;
    
    
    fprintf(fid,' %f ',Omega(i));
    for k=1:nfreq 
        w0 = Vec(1:3,k);
        wc = Vec(4:6,k);
        ws = Vec(7:9,k);
        A0_amp     = abs(w0);
        ABW_amp    = 0.5*abs(wc-sqrt(-1)*ws);
        AFW_amp    = 0.5*abs(wc+sqrt(-1)*ws);
        A0_arg     = angle(w0);
        ABW_arg    = angle(wc-sqrt(-1)*ws);
        AFW_arg    = angle(wc+sqrt(-1)*ws);
        g_amp      = abs(Vec(10:15,k));
        g_arg      = angle(Vec(10:15,k));
        
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
fclose(fid);



% Pick last 11 vectors
Freq = Freq(:,1:nfreq);

tmp = [Omega',Freq];
save('freq_struct.dat','tmp','-ascii')

stop

% Plot
figure
plot(Omega,Freq,'-o')
legend('Tower fore-aft','Tower side-side','Sym. edge/DT','BW flapwise','Sym. flap','FW flapwise','BW edgewise','FW edgewise','BW flapwise 2','Sym. flap 2','FW flapwise 2')
xlabel('Rotor speed [rad/s]')
ylabel('Frequency [Hz]')
axis([0,1.5,0,4.5])




% ------------------------------------------------ %
%% Save results for delivery

m = load('modes.dat');

% Normalize
for i=1:nfreq
    id = 1:2:29;
    m(:,30*i+1-30+id) = m(:,30*i+1-30+id)./max(max(m(:,30*i+1-30+id)));
end
save('modes_norm.dat','m','-ascii');


modeno = 5;

col1 = modeno*30-28;


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
title('Nacelle roll/torsion')
legend('Roll','Tors.')
ylabel('Amplitude')
yyaxis right
id = id+1;
plot(o,m(:,col1-1+id),'-o')
plot(o,m(:,col1-1+id+2),'-x')
ylabel('Phase')


% ------------------------------------------------ %
%% Animation

modeno = 5;
omegano = 5;

% Mode shapes
N = 50; % Number of points per mode shape
r = linspace(0,50,N); % Blades
% Blade bending
phib = BendingModeFun(r,1);
% Blade torsion
phit = TorsionalModeFun(r,1);

phi = [phib;phib;phit];


% Amplitudes
ModeAmp = abs(ModeVec(:,modeno,omegano));
% Phases
ModePhase = angle(ModeVec(:,modeno,omegano));

% Amplitude scale
AmpScale = 20;

[xb,yb,zb] = BladeGen(r);

% A single blade
Amps = ModeAmp(1:3).*phi;
Phase= ModePhase(1:3);

fig=figure;
t = 0;
dt = 0.1;
c = 2; % Chord length
while t < 1000
    t = t+dt;
    f = AmpScale*Amps(1,:).*cos(Omega(omegano)*t+Phase(1,:));
    e = AmpScale*Amps(2,:).*cos(Omega(omegano)*t+Phase(2,:));
    rot = AmpScale*Amps(3,:).*cos(Omega(omegano)*t+Phase(3,:));

    
    
    
    figure(fig);
    %plot3(r,x+0.5*rot*c,y+c/2*cos(rot));hold on
    %plot3(r,x-0.5*rot*c,y-c/2*cos(rot));hold off

    surf(xb+e,yb+f,zb,zeros(size(zb)))
        axis([-max(AmpScale*Amps(1,:))-1,max(AmpScale*Amps(1,:))+1,...
        -max(AmpScale*Amps(2,:))-2*c,max(AmpScale*Amps(2,:))+2*c,min(r),max(r),])
    
    view([30,60])
    
    pause(0.01)
    
end


