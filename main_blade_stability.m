
% DTU course main file

% Team magenta: Jan-Tore & Joey

clc
clear
close all

% ------------------------------------------------ %
% Define parameters
par = Parameters();
Wsp = 5:25;
r   = 0:par.R;
Omega = 1.5;

% ------------------------------------------------ %
% Eigen value analysis as function of omega

for i = 1:length(Wsp)

    % Aerodynamics
    AeroMat = AerodynMatrix(par,Omega,Wsp(i));
    
    
    
    % Find matrices
    [Mb,Kb,Gb,Mt,Gt,Kt]=DiagonalStructuralMatrices(par,Omega);
    ndof = size(Mb,1);

    % Define matrix A
    Ab = [zeros(ndof),  eye(ndof); -Mb\(Kb+AeroMat.Kab) , -Mb\(Gb+AeroMat.Cab)];      
   
    % Eigenvalue analysis
    [Vb,Db]=eig(Ab);
    
    % Eigen values
    lambda = diag(Db);
    [tmp,I] = sort(imag(lambda));
    lambda=lambda(I);

    id = tmp>0;
    
    lambda = lambda(id);
    Freq(i,:)=tmp(id)/2/pi; % Pick positive frequencies
    Damp(i,:)= -real(lambda)./abs(lambda);
    
    % Eigen vectors
    Vb             = Vb(:,I); % Sort vectors according to ascending frequencies
    Vec            = Vb(1:ndof,(ndof+1):(2*ndof)); % Pick vectors
    ModeVec(:,:,i) = Vec;
    
    
end



% ------------------------------------------------ %
% Save results for delivery

% Save frequencies and damping - ascending
tmp = [Freq,Damp];
tmp = [Wsp',tmp(:,[1,4,2,5,3,6])];
save('freq.dat','tmp','-ascii')

figure; hold on
plot(Wsp,Freq(:,1),'-ok')
plot(Wsp,Freq(:,2),'-xk')
plot(Wsp,Freq(:,3),'-dk')
xlabel('Wind speed [m/s]')
ylabel('Frequency [Hz]')
legend('Flapwise','Edgewise','Torsion')
yyaxis right
ylabel('Damping ratio [-]')
plot(Wsp,Damp(:,1),'-o')
plot(Wsp,Damp(:,2),'-x')
plot(Wsp,Damp(:,3),'-d')


% Define blade mode shapes
phib = BendingModeFun(r,1);
phit = TorsionalModeFun(r,1);
phi  = permute([phib;phib;phit],[1,3,2]);


% Plot mode shapes
figure
plot(r,permute(phi,[3,1,2]))
xlabel('omega [rad/s]')


% Multiply with deflefction
Amp_stand = ModeVec(:,:,4).*phi;
% Normalize
Amp_stand_norm = Amp_stand./max(max(abs(Amp_stand),[],3),[],1);
% Phase
Arg_stand = angle(Amp_stand_norm);
% Save to file
Amp=abs(permute(Amp_stand_norm,[3,1,2]));
Arg=permute(Arg_stand,[3,1,2]);
tmp = [Amp,Arg];
tmp = tmp(:,[1,4,2,5,3,6],:);
tmp = reshape(tmp,length(r),size(tmp,2)*size(tmp,3));
tmp = [r',tmp];
save('modes_8ms.dat','tmp','-ascii');

% Plot results
figure('name','8ms');hold on
subplot(3,2,1)
plot(r,Amp(:,:,1),'-o')
legend('Flap.','Edge.','Tors.')
subplot(3,2,2)
plot(r,Arg(:,:,1),'-o')
subplot(3,2,3)
plot(r,Amp(:,:,2),'-o')
subplot(3,2,4)
plot(r,Arg(:,:,2),'-o')
subplot(3,2,5)
plot(r,Amp(:,:,3),'-o')
xlabel('r [m]')
ylabel('Amplitude')
subplot(3,2,6)
plot(r,Arg(:,:,3),'-o')
xlabel('r [m]')
ylabel('Phase')

% Multiply with deflefction
Amp_stand = ModeVec(:,:,12).*phi;
% Normalize
Amp_stand_norm = Amp_stand./max(max(abs(Amp_stand),[],3),[],1);
% Phase
Arg_stand = angle(Amp_stand_norm);
% Save to file
Amp=abs(permute(Amp_stand_norm,[3,1,2]));
Arg=permute(Arg_stand,[3,1,2]);
tmp = [Amp,Arg];
tmp = tmp(:,[1,4,2,5,3,6],:);
tmp = reshape(tmp,length(r),size(tmp,2)*size(tmp,3));
tmp = [r',tmp];
save('modes_16ms.dat','tmp','-ascii');


% Plot results
figure('name','16ms');hold on
subplot(3,2,1)
plot(r,Amp(:,:,1),'-o')
legend('Flap.','Edge.','Tors.')
subplot(3,2,2)
plot(r,Arg(:,:,1),'-o')
subplot(3,2,3)
plot(r,Amp(:,:,2),'-o')
subplot(3,2,4)
plot(r,Arg(:,:,2),'-o')
subplot(3,2,5)
plot(r,Amp(:,:,3),'-o')
xlabel('r [m]')
ylabel('Amplitude')
subplot(3,2,6)
plot(r,Arg(:,:,3),'-o')
xlabel('r [m]')
ylabel('Phase')

