%***************************************************************************

% MATLAB version of the Johns [DAO, 1988] 1-D Baroclinic Instability Model 

% It employs the idealized buoyancy frequency and mean velocity
%     vertical profiles according to Gill et al. [DRS, 1974].

%                             by I.C.A.S. 

%***************************************************************************

close all
format short

% The Gill et al. [1974] parameters 

lat=input('Enter Latitude of 20 or 30 N: ');

f0=2*7.29e-5*sin(lat*pi/180);               % Coriolis parameter
beta=2*7.29e-5*cos(lat*pi/180)/6371e3;      % Rossby parameter

H=4500;  % total depth
dz=25;   % depth increment

z=-(-2*dz:dz:(H+2*dz))';
dN=900;                       % e-folding depth scale for the buoyancy frequency profile
N2=f0*f0*1e4*exp(z/dN);       % buoyancy frequency SQUARED profile
N=sqrt(N2);                   % buoyancy frequency profile
Nz=-gradient(N,dz);           % buoyancy frequency vertical gradient profile


% Gill et al. [1974] Velocity profiles are next
% There are 5 study cases to explore. Profiles [1],[2] and [3] are
% unstable for westward currents. Profile [4] is stable. Profile [5]
% is unstable for eastward currents only (Umax > 0).  

kase=input('Enter the Gill et al. [1974] study case (1 to 5): ');

Umax=-0.05;    % ***change Umax sign if desired***

if kase==1,
   % set U1=0  and U0=Umax  for study case [1] 
   U0=Umax; U1=0; dU=100;

else

 if kase==2,
   % set zmax=0 and dU=100     for study case [2]    
   zmax=0; dU=100;

elseif kase==3,
   % set zmax=-100 and dU=100  for study case [3]    
   zmax=-100; dU=100;

elseif kase==4,
   % set zmax=0 and dU=400     for study case [4]    
   zmax=0; dU=400;

else
   % set zmax=-100 and dU=400  for study case [5]    
   zmax=-100; dU=400; 
 end

% calculate U0 and U1 for case neq 1

   U0=Umax*exp(-zmax/dN)*dN/(dN-dU);
   U1=Umax*exp(-zmax/dU)*dU/(dN-dU);

end

U=U0*exp(z/dN)- U1*exp(z/dU);  % velocity profile
Uz=-gradient(U,dz);            % velocity vertical gradient profile


% compute the potential vorticity meridional gradient porfile Qy

str=-gradient(f0*f0*Uz./N2,dz);
Qy= beta - str;

% eliminate unreal data

nz=length(z);
z=z(3:nz-2);
Qy=Qy(3:nz-2);
N=N(3:nz-2);
N2=N2(3:nz-2);
Nz=Nz(3:nz-2);
U=U(3:nz-2);
Uz=Uz(3:nz-2);

% plotting N and Nz

figure(1)
subplot(121)
plot(N,z,'r','linewidth',2)
title('Buoyancy frequency profile')
xlabel('N in rad s^{-1}')
ylabel('depth in meters')

subplot(122)
plot(Nz,z,'b','linewidth',2)
title('Buoyancy frequency vertical gradient profile')
xlabel('dN/dz in rad (ms)^{-1}')
ylabel('depth in meters')

% plotting U and Uz

figure(2)
subplot(121)
plot(U,z,'r','linewidth',2)
hold on
plot([0 0],[0 -H],'k')
hold off 
title('Background Velocity Vertical profile')
xlabel('U in m s^{-1}')
ylabel('depth in meters')

subplot(122)
plot(Uz,z,'b','linewidth',2)
hold on
plot([0 0],[0 -H],'k')
hold off
title('Background Velocity Vertical Gradient profile')
xlabel('dU/dz in s^{-1}')
ylabel('depth in meters')


% plotting Qy

figure(3)
plot(Qy,z,'g','linewidth',2) 
hold on
plot([0 0],[0 -H],'k')
axis([1e-11*floor(1e11*min(Qy)) 3e-11 -H 0])
title('Potential Vorticity Meridional Gradient')
xlabel('Qy in (m s)^{-1}')
ylabel('depth in meters')
hold off

% plotting for presentation

figure(10)
subplot(121)
plot(N,z,'r','linewidth',2)
title('Buoyancy frequency profile')
xlabel('N in rad s^{-1}')
ylabel('depth in meters')

subplot(122)
plot(U,z,'b','linewidth',2)
hold on
plot([0 0],[0 -H],'k')
hold off 
title('Background Velocity Vertical profile')
xlabel('U in m s^{-1}')
ylabel('depth in meters')

figure(11)

subplot(121)
plot(Uz,z,'c','linewidth',2)
hold on
plot([0 0],[0 -H],'k')
hold off
title('Background Velocity Vertical Gradient profile')
xlabel('dU/dz in s^{-1}')
ylabel('depth in meters')

subplot(122)
plot(Qy,z,'g','linewidth',2) 
hold on
plot([0 0],[0 -H],'k')
axis([1e-11*floor(1e11*min(Qy)) 3e-11 -H 0])
title('Potential Vorticity Meridional Gradient')
xlabel('dQ/dy in (m s)^{-1}')
ylabel('depth in meters')
hold off

pause

%----------------------------------------------------------------------

% The Eigenvalue Problem

% Based on the original Johns [1988] fortran 77 code 

% Consists of solving complex eigenvalues (c) an eigenvectors (P)

% The QG vertical equations AND boundary conditions are discretized 
% using FULL CENTERED FINITE DIFFERENCES

% The discretized system has the form (A-Bc)P=0 of nz linear equations
% with the unknowns P(n), n=1,2..nz.

% A and B are tridiagonal matrices

%----------------------------------------------------------------------

nz=length(z);   % number of equations/levels

% The disturbance amplitude channel is particular to John's model,
% not Gill et al.'s. Simply use L=Inf to zero the channel part out of
% the problem

L=input('Enter the Disturbance Amplitude Channel in km: ');

L=1e3*L;

% the bottom slope parameter is scaled by the slope of the deep isopycnals
% Gill et al. [1974 considered the scaled s=-10,0,10
% this |s=10| value corresponds to a nonscaled value of about s=1.12e-3

s=input('Enter Bottom Slope: ');
display(' bottom slope scaled value is:')
s=s*(N2(nz)/Uz(nz)/f0)  % scaling bottom slope


k=0.05:-0.001:0.001; % wavenumbers in km^(-1)
k=k*1e-3;
nk=length(k);

% initialize the phase speed, growth rate and most unstable mode matrices
cr=zeros(nk,1);
sig=cr;
P=zeros(nz,nk);
teste=zeros(nk,1);

for n=1:nk,                 % begin wavenumber loop

    n

% set up matrix coefficients

    L1= f0*f0./N2.*( (1/dz/dz) - (1/dz)*Nz./N );
    L2= 2*f0*f0./N2.*(1/dz/dz) + (k(n)*k(n) + pi*pi/L/L);
    L3= f0*f0./N2.*( (1/dz/dz) + (1/dz)*Nz./N );

% build tridiagonal matrices A and B

A=diag(Qy-L2.*U)+diag(L3(1:(nz-1)).*U(1:(nz-1)),1)+diag(L1(2:nz).*U(2:nz),-1);

B=diag(-L2)+diag(L3(1:(nz-1)),1)+diag(L1(2:nz),-1);

% fix extreme points with Boundary Conditions

A(1,1)=A(1,1)+L1(1)*Uz(1)*2*dz;        % top BC
A(1,2)=A(1,2)+L1(1)*U(1);

A(nz,nz)=A(nz,nz)-L3(nz)*(1-s)*Uz(nz)*2*dz;  % bottom BC
A(nz,nz-1)=A(nz,nz-1)+L3(nz)*U(nz);

B(1,2)=B(1,2)+L1(1);                  % top BC
B(nz,nz-1)=B(nz,nz-1) + L3(nz);       % bottom BC

% obtain the eigenvalue matrix

C=B\A;

% calculate eigenvalues and eigenvectors for k(n)

[F,lamb]=eig(C);

% save most unstable mode and growth rates corresponding to k

[ci,jmax]=max(diag(imag(lamb)));

        if ci ==0,
sig(n)=0;
cr(n)=NaN;
P(:,n)=NaN*ones(nz,1);

else
sig(n)=k(n)*ci*86400;              % growth rate in days^(-1)
cr(n)=real(lamb(jmax,jmax))*100;   % phase speed in cm/s
P(:,n)=F(:,jmax); 
teste(n)=k(n);                 % vertical mode
        end

end                         % end wavenumber loop                                      


% compute amplitude and phase of the most unstable modes

Pamp=abs(P);                         % amplitude
Pphase=180/pi*angle(P);              % phase in degrees

% find and normalize most unstable mode

[sigmax, imax]=max(sig);

Pmax=Pamp(:,imax); 

Pmax=Pmax/max(abs(Pmax));  % normalization by maximum value
%Pmax=Pmax/(norm(Pmax));    % orthonormalization 

Pphmax=Pphase(:,imax);


% plot phase speeds and growth rates

figure(4)

subplot(211)

plot(k*1000,-cr,'b','linewidth',2)
axis([0 0.05 0 6])
title('Phase Speeds')
ylabel('phase speeds in cm^{-1}')
xlabel('wavenumber in km^{-1}') 

subplot(212)

plot(k*1000,sig,'r','linewidth',2)
axis([0 0.05 0 0.016])
title('Phase Speeds')
ylabel('growth rates in days^{-1}')
xlabel('wavenumber in km^{-1}') 

% plot amplitude and phase of the most unstable modes

figure(5)
subplot(121)
plot(Pmax,z,'g','linewidth',2)
hold on
plot([0 0],[0 -H],'k')
hold off 
title('Most Unstable Mode Amplitude')
xlabel('P mode amplitude')
ylabel('depth in meters')

subplot(122)
plot(Pphmax-Pphmax(nz),z,'m','linewidth',2)  % 
%axis([-90 30 -H 0])
hold on
plot([0 0],[0 -H],'k')
hold off
title('Phase relative to the bottom value')
xlabel('Phase in degrees')
ylabel('depth in meters')

