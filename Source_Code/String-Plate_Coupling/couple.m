%-------------------------------------------------------------------------%
% FDTD Model for string plate connection. The connection here is assumed
% to be rigid with and without loss.
%Author:Mangesh Chandrakant Sonawane
%MSc Acoustic and Music Technology dissertation project
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
                                % Flags
%-------------------------------------------------------------------------%
 
 bctype_p = 1; %1: simply supported, 2: clamped
 bctype_S = 1; %1: simply supported, 2: clamped
 outtype=1;    % output type: 1: displacement, 2: velocity
 losstype = 2; % loss type: 0:no loss 1: independant, 2: dependant
 intype = 1;    % type of input: 1: struck, 2: plucked
 energyAn=1;

%-------------------------------------------------------------------------%
                        % Simulation Parameters
%-------------------------------------------------------------------------%

 % simulation
 Tf = 2; % duration in seconds
 % readout position as percentage
 readout_position = [.85 .85; .95 .95]; 
% sample rate (Hz)
 SR = 44.1e3; 
 k = 1/SR;
%-------------------------------------------------------------------------%
                            % Physical Parameters
%-------------------------------------------------------------------------%
% Plate

% Young's modulus
 E_Plate = 2e11;
% density (kg/m^3)
 rho_Plate =7850;
% Poisson Ratio (< .5)
 nu_Plate=.3;
% thickness (m)
 H_Plate = .005;
% plate length X axis (m)
 Lx_Plate = .5;
% plate length X axis (m)
 Ly_Plate = .4;
% loss [freq.(Hz), T60;...]
 loss_Plate = [100, 2; 1000, 1.5];
% coordinates of string coupling [Xs1,Ys1; Xs2,Ys2];
 ctr_Plate = [.35 .35; .71 .78];

% String

%Frequency for note
 f0_String = 440; %A6
% string radius (m)
 r_String = 3e-4;
% string lengths 
 L_String = 1;
% Young's modulus 
 E_String = 2e11;%6e11;
% density (kg/m^3)
 rho_String = 7850;%1350;
% loss [freq1(Hz), T60_1; freq_2, T60_2] 
 loss_String = [100, 6; 900 5];
% I/O string
 xi_String = 0.75;  % coordinate of excitation (normalised, 0-1)
% peak amplitude of excitation (N)
 famp_String = 1;
% duration of excitation (s)
 dur_String = 0.001; 
% start time of excitation (s) 
 st_exc_st = 0.04;
 % Tension in Newtons using frequency of particular note
 st_T = (((2*f0_String.*r_String).*L_String).^2)*pi*rho_String; 
 %st_T=99.279e1; %Alternate for frequency
 
%-------------------------------------------------------------------------%
                        % Plate Derived Parameters
%-------------------------------------------------------------------------%

 % Coefficients
 D_Plate = (E_Plate*(H_Plate)^3)/(12*(1-(nu_Plate^2)));
 kappa_Plate = sqrt(D_Plate / (rho_Plate*H_Plate) );
 hmin_Plate = 2*sqrt(k*kappa_Plate);
 Nx_Plate = floor((Lx_Plate)/hmin_Plate);
 Ny_Plate = floor((Ly_Plate)/hmin_Plate);
 %Grid Spacing
 h_Plate = sqrt(Lx_Plate*Ly_Plate)/sqrt(Nx_Plate*Ny_Plate);
 mu_Plate = (kappa_Plate * k)/(h_Plate^2);
 % % number of time steps
 Nf = floor(SR*Tf);
 Nx_Plate = Nx_Plate +1;
 Ny_Plate = Ny_Plate +1;
% total grid 
 tg = Nx_Plate*Ny_Plate; 
 
%-------------------------------------------------------------------------%
                        % String Derived Parameters
%-------------------------------------------------------------------------%

% Cross-sectional area
 A_String = pi*r_String^2;
% Moment of intertia
 I_String = 0.25*pi*r_String^4;
 % wave speed
 c_String = sqrt(st_T/(rho_String*A_String));
 % stiffness constant
 K_String = sqrt(E_String*I_String/(rho_String*A_String));
 % minimal grid spacing for stability
 hmin_String = sqrt(0.5* (c_String.^2*k^2+sqrt(c_String.^4*...
     k^4+16*K_String.^2.*k.^2)) );
% number of grid points to update
 N_String = floor(L_String/hmin_String); 
% actual grid spacing used
 h_String = L_String/N_String;
% Courant number (?)
 lambda_String = c_String*k/h_String;
% numerical stiffness constant (?) 
 mu_String = K_String*k/h_String^2;
%Grid points update
 N_String = N_String +1;

 
%-------------------------------------------------------------------------%
                                  % Loss
%-------------------------------------------------------------------------%

 if losstype==0
 sigma0_String=0;
 sigma1_string=0;
 sigma0_Plate=0; 
 sigma1_Plate=0;
 end

 if losstype==1
 % frequency independant loss
 sigma0_String = 6*log(10)/loss_String(1,2);
 sigma1_string = 0;
 sigma0_Plate = 6*log(10)/loss_Plate(1,2);
 sigma1_Plate = 0;
 end

 if losstype == 2
    
 
     % String
 coeff1_String = (-c_String.^2 + sqrt(c_String.^4 + 4*K_String.^2.*...
     (2*pi*loss_String(1,1))^2))./(2*K_String.^2);
coeff2_String = (-c_String.^2 + sqrt(c_String.^4 +4*K_String.^2.*...
    (2*pi*loss_String(2,1))^2))./(2*K_String.^2);

 sigma0_String = 6*log(10)*(-coeff2_String/loss_String(1,2)...
     + coeff1_String/loss_String(2,2))./(coeff1_String -coeff2_String);
 sigma1_string = 6*log(10)*(1/loss_String(1,2)...
     - 1/loss_String(2,2))./(coeff1_String-coeff2_String);
     
% Plate
 coeff1_Plate = 2*kappa_Plate*(2*pi*loss_Plate(1,1))/(2*kappa_Plate.^2); 
 coeff2_Plate = 2*kappa_Plate*(2*pi*loss_Plate(2,1))/(2*kappa_Plate.^2);
 sigma0_Plate = 6*log(10)*(-coeff2_Plate/loss_Plate(1,2) + coeff1_Plate/...
     loss_Plate(2,2))./(coeff1_Plate -coeff2_Plate);
 sigma1_Plate = 6*log(10)*(1/loss_Plate(1,2) - 1/loss_Plate(2,2))...
     ./(coeff1_Plate-coeff2_Plate);




 end 
%-------------------------------------------------------------------------%
                            % Read In/Out
%-------------------------------------------------------------------------%
 lo = readout_position .*[ Nx_Plate Ny_Plate ];
 lo = [floor(sub2ind([Nx_Plate Ny_Plate],lo(1), lo(3))),...
     floor(sub2ind([Nx_Plate Ny_Plate],lo(2), lo(4)))];
%-------------------------------------------------------------------------%
                         % Create Force Signal
%-------------------------------------------------------------------------%

 f_String = zeros(Nf,1);
 durint = floor(dur_String*SR);
 exc_st_int = (floor(st_exc_st*SR))+1; 
 durf = exc_st_int:exc_st_int+durint -1;
 fcoeff_string = (k^2/(h_String*rho_String*A_String));
 f_String(durf) = famp_String*0.5*(1-cos((2/intype)*pi.*(durf...
     /durint)))*fcoeff_string;

%-------------------------------------------------------------------------%
                          % Plate Coefficient 
%-------------------------------------------------------------------------%
%Laplacian
laplacian = triangle(Ny_Plate ,Nx_Plate ,1,bctype_p); 
%Biharmonic
biharmonic = triangle(Ny_Plate ,Nx_Plate ,2,bctype_p); 
I = speye(tg);
I([1 end],:) = 0;
out1 = I;
pl_mI = out1;
A_coeff_Plate =(1/(1+k*sigma0_Plate))*speye(tg);
B_coeff_Plate =(-(mu_Plate^2)*biharmonic + 2*pl_mI + (2*k*...
    sigma1_Plate/h_Plate^2)*laplacian) * A_coeff_Plate;
C_coeff_Plate = (-(1-sigma0_Plate*k)*pl_mI - (2*k*sigma1_Plate/...
    h_Plate^2)*laplacian) * A_coeff_Plate; 
%-------------------------------------------------------------------------%
                        % String Coefficients
%-------------------------------------------------------------------------%

% scheme coefficients
ls=N_String;
outMx = spdiags([-ones(ls,1),ones(ls,1)],0:1,speye(ls));
Dn=outMx;
%Second order operator
Dnn = Delta(N_String,'xx');
%Fourth order operator
Dnnnn = Delta(N_String,'xxxx',bctype_S);
Dnnnn(2,1) = -2;
st_mI = speye(N_String);
% String matrix 
s_A = (1/(1+k*sigma0_String));
s_B = ((lambda_String^2)*Dnn - (mu_String^2)*Dnnnn + (2*st_mI) ...
     + ((2*sigma1_string*k/h_String^2)*Dnn)) * s_A ;
s_C = -((1-sigma0_String*k)*st_mI + ((2*sigma1_string*k/h_String^2)...
     *Dnn) ) * s_A;
%coupling string boundary
 s_B([1 end],:) = 0;
 s_C([1 end],:) = 0;
 s_B(1,:) = ((lambda_String^2)*Dn(1,:) - (mu_String^2)*Dnn(2,:)...
     + (2*st_mI(1,:))) * s_A;
 s_C(1,:) = -((1-sigma0_String*k)*st_mI(1,:));

%-------------------------------------------------------------------------%
                  % Matrices update and Couple Matrices
%-------------------------------------------------------------------------%

 %Total Number of Grid points
 Nt = tg+ N_String;
% coupling points
 w0p = sub2ind([Nx_Plate, Ny_Plate], floor(ctr_Plate(1,1)...
     *Nx_Plate),floor(ctr_Plate(1,2)*Ny_Plate));
% Mass Ratio Coefficients
 M_Plate = 1/((rho_Plate*H_Plate*h_Plate^2)*(1+k*sigma0_Plate));
 M_String = 1/((rho_String*A_String*h_String)*(1+k*sigma0_String));
% Spreading operators
 J = sparse(zeros(Nt,1)); 
 J([w0p tg+1]) = [M_Plate -M_String];
 pJ = sparse(tg,1);
 pJ(w0p) = 1;
 sJ = sparse(N_String ,1);
 sJ(1) = 1;
% Mass ratio
 M_Coef1 = 1/( 1/(rho_String*A_String*h_String*(1+k*sigma0_String)) +...
 (1/(rho_Plate*H_Plate*(h_Plate^2)*(1+k*sigma0_Plate)))*(pJ'*pJ));
% Coupling Vector
 F = M_Coef1*[-pJ'*B_coeff_Plate,sJ'*s_B];
 F1 = M_Coef1*[pJ'*(2*pl_mI + (2*sigma1_Plate*k/h_Plate^2)*laplacian)...
     /(1+k*sigma0_Plate),...
-2*sJ'*st_mI/(1+k*sigma0_String)];
%Sparsed Matrix as coefficient for main calculation
 B=blkdiag(B_coeff_Plate ,s_B) + J*F;
 C=blkdiag(C_coeff_Plate ,s_C) + J*F1;

%-------------------------------------------------------------------------%
                             % Initialise I/O
%-------------------------------------------------------------------------%
% initialise output
u =zeros(tg,1); 
w=zeros(N_String ,1);
% Joined vectors
 uw = [u;w];
 uw1 = uw;
 uw2 = uw;
% input
 fvect = [u;w];
% Index of excitation
 li = floor(xi_String*N_String) + tg; 
% output
 y = zeros(Nf,2);
%Energy analysis initialise
if energyAn==1
  % total energy
  plate_Energy = zeros(Nf,1); 
  % kineteic energy
  plate_KE =zeros(Nf,1);  
  % potential enery
  plate_PE = zeros(Nf,1);
  % total energy
  string_Energy = zeros(Nf,1); 
  % kineteic energy
  string_KE = zeros(Nf,1); 
  % potential energy
  string_PE = zeros(Nf,1);        
  fEnergy = zeros(Nf,1);
end


 
%-------------------------------------------------------------------------%
                            % Main Calculation
%-------------------------------------------------------------------------%

 tic

 for n = 1:Nf 
 % update input forces
 fvect(li) = f_String(n);

 % main operation
 uw = B*uw1 + C*uw2 + fvect;

 % read output
 if (outtype==1)
 y(n,:) = uw(lo);

 elseif (outtype==2)
 y(n,:) = (SR*(uw(lo)-uw1(lo)));

 end 

 
 if energyAn==1 
      u = uw(1:tg); 
      u1 = uw1(1:tg);
      w = uw(tg+1:end); 
      w1 = uw1(tg+1:end);

      plate_KE(n) = 0.5*(h_Plate^2)/k^2*((u-u1)'*(u-u1));
      plate_PE(n) = 0.5*(kappa_Plate^2)/(h_Plate^2)*((laplacian*u)'...
          * (laplacian*u1));
      plate_Energy(n)=plate_KE(n)+plate_PE(n);
     % pEnergy(n) = pl_coE*[sum(((u-u1).^2)); (laplacian*u)' * (laplacian*u1)];

      string_KE(n) = 0.5*(h_String)/k^2 *((w-w1)'*(w-w1));
      string_PE(n) = 0.5*(c_String^2)/h_String*((Dn*w)' * (Dn*w1));
      string_Energy(n) = string_KE(n)+string_PE(n);
      %fEnergy(n)=pEnergy(n)+sEnergy(n);
      fEnergy(n) = (k*c_String^2/(2*h_String^2))*(w(2)-w(1)) ...
          + SR*(w(1)-w1(1)) + (kappa_Plate^2*k*.5)*B_coeff_Plate(w0p,:)...
          *u - SR*(u(w0p) - u1(w0p));
 end
   
 
 % Update
 uw2=uw1;
 uw1=uw;

 end 
 y= y/max(abs(y));
soundsc(y,SR);
 toc
%  n=1:Nf;
%  plot(n,y(n));
%  xlabel('Samples');
%  ylabel('magnitude');

% figure
%       yfft = 10*log10(abs(fft(y)));
%       plot([0:Nf-1]'/Nf*SR, -abs(yfft), 'k')
%       xlabel('freq (Hz)')
%       ylabel('|magnitude| (db)');
%       title('total system output stectra')
% 
% plot(y);
% xlabel('Sample');
% ylabel('Amplitude');
% title('String-Plate connection displacement');

