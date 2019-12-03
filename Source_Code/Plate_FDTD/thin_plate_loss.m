%-------------------------------------------------------------------------%
% A Kirchhoff Thin Plate FDTD Model with and without loss
%Author:Mangesh Chandrakant Sonawane
%MSc Acoustic and Music Technology dissertation project
%-------------------------------------------------------------------------%
clear all
close all

% boundary condition type: 1: simply supported, 2: clamped
bctype = 1;
% output type: 1: displacement, 2: velocity
outtype = 1;  
% loss type: 0: independent, 1: independent 2: Freq dependent
losstype = 1; 

%-------------------------------------------------------------------------%
                    % Parameters for simulation
%-------------------------------------------------------------------------%
% sample rate (Hz)
SR = 44100; 
% duration
Tf = 10;  
% number of time steps
Nf = floor(SR*Tf); 
% Poisson Ratios (< .5)
poissnratio = 0.5; 
%excitation (normalised)
exc = [0.45, 0.45];             
wi = 0.2; 
%initial displacement
u0 = 0;
%initial velocity
v0 = 1;    
% readout position as percentage 
rp = [.45, .65];
%-------------------------------------------------------------------------%
                        % physical parameters
%-------------------------------------------------------------------------%
% Young's modulus
E = 2e11; 
% density (kg/m^3)
rho = 7850; 
% thickness (m)
H = .005; 
% plate length (m)
L = .9; 
% plate length X axis (m)
Lx = .7;      
% plate length Y axis (m)
Ly = .7;                     
 % loss [freq.(Hz), T60;...]
loss = [500, 2; 1500, 1];  
 %  Motion Coefficients
D = (E*(H)^3)/(12*(1-(poissnratio^2)));
kappa = sqrt(D / (rho*  H) );
%Time step
k = 1/SR;   
%Spacing
hmin = 2*sqrt(k*kappa);      
N = floor(L./hmin); 
h = L./(N); 
%Nx = floor(Lx/hmin)-1
%Ny =floor(Ly/hmin)-1
%h=sqrt(Lx*Ly)/sqrt(Nx*Ny);
mu = (kappa * k)/(h^2);
%Updated value
N = N+1;
% total grid size.
ss = N*N;                   
%-------------------------------------------------------------------------%
                        % Loss coefficients
%-------------------------------------------------------------------------%
if losstype ==0
  % no loss
  sigma0 = 0;
  sigma1 = 0;
end

if losstype ==1
  % frequency independant loss 
  sigma0 = 6*log(10)/loss(1,2);
  sigma1 = 0;
end

if losstype == 2
  %frequency  dependent
  z1 = 2*kappa*(2*pi*loss(1,1))/(2*kappa.^2); 
  z2 = 2*kappa*(2*pi*loss(2,1))/(2*kappa.^2);

  sigma0 = 6*log(10)*(-z2/loss(1,2) + z1/loss(2,2))./(z1-z2);
  sigma1 = 6*log(10)*(1/loss(1,2) - 1/loss(2,2))./(z1-z2);
end
%-------------------------------------------------------------------------%
                          % Read In/Out
%-------------------------------------------------------------------------%
lo = rp*N;
lo = floor(sub2ind([N N],lo(1), lo(2)));
li = exc*N;
li = floor(sub2ind([N N],li(1), li(2)));
%-------------------------------------------------------------------------%
                        % Create Force Signal
%-------------------------------------------------------------------------%
[X,Y] = meshgrid([1:N]*h, [1:N]*h);     
% distance of points from excitation
dist = sqrt((X-(exc(1)*L)).^2 + (Y-(exc(2)*L)).^2); %NSS
ind = sign(max(-dist+(wi*0.5),0));         
rc = 0.5*ind.*(1+cos(2*pi*dist/wi));        %NSS
rc = rc(:);                                 
%-------------------------------------------------------------------------%
                              %operators
%-------------------------------------------------------------------------%
biharmonic_ope = triangle(N,N,2,bctype); % biharmonic matrix operator=2
laplacian_ope = triangle(N,N,1,bctype);  % Laplacian matrix operator=1
%-------------------------------------------------------------------------%
                        % Coefficient Matrices
%-------------------------------------------------------------------------%
A = (1/(1+k*sigma0))*speye(ss);   
B = (-(mu^2)*biharmonic_ope + (2*sigma1*k/(h^2))*laplacian_ope + 2*speye(ss)) * A;
C = (-(2*sigma1*k/(h^2))*laplacian_ope - (1-sigma0*k)*speye(ss))  * A;
%-------------------------------------------------------------------------%
                            % initialise output
%-------------------------------------------------------------------------%
u2 = u0*rc;
u1 = (u0+(k*v0))*rc;
u  = u2;
y = zeros(Nf,1);
%-------------------------------------------------------------------------%
                                 % Main Calculation
%-------------------------------------------------------------------------%

  for n = 1:Nf

    % main operation
    u = B*u1 + C*u2;
    % read output
    if (outtype==1)
      y(n,:) = u(lo);
    elseif (outtype==2)
      y(n,:) = (SR*(u(lo)-u1(lo)));
    end
     u2 = u1; u1 = u;
  end

y= y/max(abs(y));

figure
%       yfft = 10*log10(abs(fft(y)));
%       plot([0:Nf-1]'/Nf*SR,-abs(yfft), 'k')
%       xlabel('freq (Hz)')
%       ylabel('|magnitude| (db)');
%       title('Plate Frequency response')
plot(y)
      xlabel('Sample');
      ylabel('Amplitude');
      title('Plate Displacement');
soundsc(y,SR)
audiowrite('Plate_Loss.wav', y,SR)
%-------------------------------------------------------------------------%
