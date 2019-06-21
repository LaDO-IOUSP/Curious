# -*- coding: utf-8 -*-

# packages
import numpy as np
import seawater as sw

def BCinst(U,N2,z,Qy,k,lat,s,L=np.inf,ct=50):

    ''' Computes 1-d baroclinic instability of a vertical velocity profile 
    ===========================================================================
    
    Input :
    -> U:     velocity profile [m/s]
    -> N2:    squared Brunt-Vaisala frequency [1/s]
    -> z:     depth (z = 0 at surface) [m] 
    -> dQ/ds: cross-stream vorticity gradient [1/m*s]
    -> k:     wavenumber array [1/km]
    -> lat:   mean latitude [-90 ... 90 deg]
    -> s:     bottom slope [unitless]
    
    * args
    -> L:  the disturbance amplitude channel [km] (default "np.inf" - problem 
           independent of a channel)
    -> ct: cut off parameter to eliminate unreal data [m] (default: 50)
    ===========================================================================
    
    Output:
    -> sig:    instability growth rates [1/days]
    -> wl:     wavelengths [km]
    -> Pmax:   amplitude of max unstable wavenumber [unitless]
    -> Pphmax: phase of max unstable wavenumber [deg]
    -> imax:   index for the most unstable mode [unitless]
    -> cr:     real part of phase speed [m/s] 
    -> z:      modified depth array [m]

    ===========================================================================

    Python version by:
        Hélio Almeida (helio.almeida@usp.br) 
        Dante C. Napolitano (dante.napolitano@usp.br) 
        @ LaDO-IOUSP in 27/10/2016
    
    based on matlab code of (adapted from the original Johns [1988] fortran) by:
        Ilson C.A. Silveira (ilson.silveira@usp.br)
        (Available at https://github.com/LaDO-IOUSP/Curious/blob/master/Matlab/gill74.m)
        
    '''

    # compute Coriolis parameter
    f0 = sw.f(lat)
    
    # Computes Uz, N, Nz
    
    diffz=np.gradient(z)
    
    if np.unique(diffz).size!=1:
        raise ValueError('z must be equally spaced')
    Uz=np.gradient(U,diffz)
    N = np.sqrt(N2)
    Nz = np.gradient(N,diffz)
    
    
    dz = np.abs(diffz)[0]

    # eliminate unreal data
    cut=int(ct/dz)
    nz=z.shape[0]
    z=z[cut:nz-cut]
    Qy=Qy[cut:nz-cut]
    N=N[cut:nz-cut]
    N2=N2[cut:nz-cut]
    Nz=Nz[cut:nz-cut]
    U=U[cut:nz-cut]
    Uz=Uz[cut:nz-cut]


    nz=z.shape[0]   # number of equations/levels

    # Length in meters
    L=1e3*L

    # scaling bottom slope
    s1=s*(N2[nz-1]/Uz[nz-1]/f0)


    # wavenumber in meters
    k=k*1e-3
    nk=k.shape[0]

    # initialize the phase speed, growth rate and most unstable mode matrices
    cr=np.zeros([nk,1])
    sig=cr.copy()
    P=np.zeros([nz,nk])+0j


    for n in np.arange(0,nk): # begin wavenumber loop
        
        print '%i'%n
        
        # set up matrix coefficients
        L1= f0*f0/N2*((1./dz/dz)-(1./dz)*Nz/N); L1=L1.squeeze()
        L2= 2*f0*f0/N2*(1./dz/dz) + (k[n]*k[n] + pi*pi/L/L); L2=L2.squeeze()
        L3= f0*f0/N2*( (1./dz/dz) + (1./dz)*Nz/N); L3=L3.squeeze()
    
        # build tridiagonal matrices A and B
    
        A=np.diag(Qy-L2*U)+np.diag(L3[0:nz-1]*U[0:nz-1],1)+np.diag(L1[1:nz]*U[1:nz],-1)
    
        B=np.diag(-L2)+np.diag(L3[0:nz-1],1)+np.diag(L1[1:nz],-1)
    
        # fix extreme points with Boundary Conditions
    
        A[0,0]=A[0,0]+L1[0]*Uz[0]*2.*dz        # top BC
        A[0,1]=A[0,1]+L1[0]*U[0]
    
        A[nz-1,nz-1]=A[nz-1,nz-1]-L3[nz-1]*(1-s1)*Uz[nz-1]*2.*dz  # bottom BC
        A[nz-1,nz-2]=A[nz-1,nz-2]+L3[nz-1]*U[nz-1]
    
        B[0,1]=B[0,1]+L1[0]                        # top BC
        B[nz-1,nz-2]=B[nz-1,nz-2] + L3[nz-1]       # bottom BC
    
        # obtain the eigenvalue matrix
    
        C=np.linalg.solve(B,A)
        # calculate eigenvalues and eigenvectors for k(n)
    
        lamb,F=np.linalg.eig(C)
    
        # save most unstable mode and growth rates corresponding to k
    
        ci=np.max(np.imag(lamb))
        jmax=np.argwhere(np.imag(lamb)==ci)[0]
    
        if ci ==0:
            sig[n,0]=0
            cr[n,0]=np.nan
            P[:,n]=np.nan*np.ones([nz])+np.nan*np.ones([nz])*1j
        else:
            sig[n,0]=k[n]*ci*86400.              # growth rate in days^(-1)
            cr[n,0]=np.real(lamb[jmax])*100      # phase speed in cm/s
            P[:,n]=F[:,jmax].squeeze()           # vertical mode
    # end wavenumber loop                                      
    
    
    # compute amplitude and phase of the most unstable modes
    
    Pamp=np.abs(P)                         # amplitude
    Pphase=180/np.pi*np.angle(P)              # phase in degrees
    
    # find and normalize most unstable mode
    
    sigmax=np.max(sig)
    imax=np.argwhere(sig==sigmax)[0][0]
    
    Pmax=Pamp[:,imax]
    
    Pmax=Pmax/np.max(np.abs(Pmax))  # normalization by maximum value
    #Pmax=Pmax/(norm(Pmax));        # orthonormalization 
    
    Pphmax=Pphase[:,imax]
    
    # obtain wavelenghts in km
    wl = np.array(0.002*np.pi/k)
    
    k = 2*np.pi/wl       
    
    if np.all(sig == 0):
        print 'No instability found in profile'
    else:    
        print 'Most unstable wavelength: %.3f' %(wl[imax])  
        
    return sig,wl,Pmax,Pphmax,imax,cr,z
