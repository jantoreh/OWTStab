
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


% ------------------------------------------------ %
% Eigen value analysis as function of omega

for i = 1:length(Omega)

    % Find matrices
    [Mb,Kb,Gb,Mt,Gt,Kt]=DiagonalStructuralMatrices(par,Omega(i));
    ndof = size(Mb,1);

    % Define matrix A
    Ab = [zeros(ndof),  eye(ndof); -inv(Mb)*Kb , -inv(Mb)*Gb];      
   
    % Eigenvalue analysis
    [Vb,Db]=eig(Ab);
    
    % Eigen values
    [tmp,I] = sort(diag(imag(Db))/2/pi); % Sort for ascending frequencies
    Freq(i,:)=tmp(tmp>0); % Pick positive frequencies
    
    % Eigen vectors
    Vb = Vb(:,I); % Sort vectors according to ascending frequencies
    Vec(:,:,i) = Vb(1:ndof,(ndof+1):(2*ndof)); % Pick vectors
    ModeAmp(:,:,i) = (Vec(:,:,i));
    ArgVec(i,:) = zeros(1,ndof);
    
end

% ------------------------------------------------ %
% Save results for delivery

% Save frequencies - ascending
tmp = [Omega',Freq];
save('freq.dat','tmp','-ascii')

% Define blade mode shapes
phib = BendingModeFun(r,1);
phit = TorsionalModeFun(r,1);
phi  = permute([phib;phib;phit],[1,3,2]);


% Plot mode shapes
figure
plot(r,permute(phi,[3,1,2]))
xlabel('omega [rad/s]')


% Multiply with deflefction
Amp_stand = ModeAmp(:,:,1).*phi;

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
save('mode_standstill.dat','tmp','-ascii');


% Multiply with deflefction
Amp_stand = ModeAmp(:,:,end).*phi;

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
save('mode_rated.dat','tmp','-ascii');






