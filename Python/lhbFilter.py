# -*- coding: utf-8 -*-
import numpy as np
def fourfilt(data,delt,T,ftype='band'):
    '''Computes data fourrier filter working as lowpass, highpass or bandpass
    ===========================================================================

    Input :
        
        -> data: the set of all data before filtering  (array Nx1 or 1xN)
        
        -> delt: sampling interval  (float)
        
        -> T: array of periods to filter data 
              T size must be 1 in case of low and highpass and
              T size 2 in case of bandpass filter(i.e T=np.array([Tmin,Tmax]))
        
        -> ftype: type of filtering to do : low, high or band (string);
      
        
        CAUTION: ftype default is bandpass
        
    ===========================================================================

    Output:
        
        -> filtdata: filtered data (same as input data (Nx1))
        
    ===========================================================================

    Python version by:
        Hélio Almeida (helio.almeida@usp.br) 
        @ LaDO-IOUSP in 04/10/2016
    
    based on matlab code of:
        Jeff List (jlist@usgs.gov) in (12/4/96)
        Rich Signell (rsignell@usgs.gov) in (1/8/97)
    '''
    
    try:
        s1,s2=data.shape
    except:
        s1=data.shape
        s2=1
    if (s1!=1)&(s2!=1):
        raise ValueError('fourfilt cannot handle matrices')

    ft=ftype.lower()

    if (ft!='low')&(ft!='high')&(ft!='band'):
        raise ValueError('Wrong ftype chosen, please choose low, high or band!')
    
    if ft=='low':
        tmin = T[0]
        tmax = (data.size*delt)+1
    elif ft=='high':
        tmin=(2*delt)-1
        tmax = T[0]
    else:
        tmin=np.min(T)
        tmax=np.max(T)
        
    npts=data.size
    
    nby2=npts/2.
    tfund=npts*delt
    ffund=1/tfund
    
    #  remove the mean from data:
    
    datamean=np.nanmean(data)
    data-=datamean
    
    # fourier transform data:
    
    coeffs=np.fft.fft(data)
    
    #  filter coefficients:
    
    f=ffund
    for i in np.arange(1,int(nby2)+1):
        t=1.0/f
        if (t > tmax) | (t < tmin):
            coeffs[i]=0.0*coeffs[i]
        f=f+ffund
    
    
    
    #  calculate the remaining coefficients:
    
    for i in np.arange(1,int(nby2)):
        coeffs[npts-i]=np.conj(coeffs[i])
    
    
    #  backtransform data and take real part:
    
    backcoeffs=np.fft.ifft(coeffs)
    filtdat=np.real(backcoeffs)
    
    # add back the mean:
    filtdat=filtdat+datamean
    
    return filtdat