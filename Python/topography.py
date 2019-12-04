# -*- coding: UTF-8 -*-
import netCDF4 as nc
import numpy as np
import scipy.interpolate as scint

def geb_bat(filename):

    data = nc.Dataset(filename)
    LON = data.variables[u'lon'][:]
    LAT = data.variables[u'lat'][:]
    DEPTH = data.variables[u'elevation'][:]
            
    DEPTH[DEPTH>0] = 0
    DEPTH = -1.*DEPTH
     
    LONb,LATb = np.meshgrid(LON,LAT) 
    
    return LONb,LATb,DEPTH
    
# Smoothing
def smoo(self,nw):
    import scipy.signal as sg
    win = np.blackman(nw)
    win=win[None,:]*win[:,None]
    peso=sg.fftconvolve(np.ones(self.shape),win,mode='same')
    dataout=sg.fftconvolve(self,win,mode='same')/peso
    
    return dataout 
    
def geb_bat_smoothed(filename,lev=50):   

    LONb,LATb,DEPTH = geb_bat(filename)
    batsmoo=smoo(DEPTH,lev)
    
    return LONb,LATb, batsmoo
    
#######################################################
# extract bat from radial
def extract_bat(lonctd,latctd,angle,PVEL,DVEL,VMDR,distmax):
    
    #A DEFINIR LIMITES DO MAPA
    lnu,lnd = lonctd.max()+1,lonctd.min()-1
    ltu,ltd = latctd.max()+1,latctd.min()-1
    
    #A LER DADOS DE BATIMETRIA GEBCO
    
    bLON,bLAT,BAT=geb_bat()
    
    #ENCONTRANDO PROFUNDIDADE PARA CADA PONTO DA RADIAL
    radctd = np.array([])
    for lonrad,latrad in zip(lonctd,latctd):
        l1 = near(bLON[0,:],lonrad,1)
        l2 = near(bLAT[:,0],latrad,1)
        radctd = np.append(radctd,BAT[np.argwhere(bLAT[:,0]==l2),
                        np.argwhere(bLON[0,:]==l1)])
    
    #ENCONTRANDO PROFUNDIDADE PARA CADA PONTO DA RADIAL
    lonmed = lonctd[:-1]+np.diff(lonctd)/2
    latmed = latctd[:-1]+np.diff(latctd)/2 
    
    radmed = np.array([])
    for lonrad,latrad in zip(lonmed,latmed):
        l1 = near(bLON[0,:],lonrad,1)
        l2 = near(bLAT[:,0],latrad,1)
        radmed = np.append(radmed,BAT[np.argwhere(bLAT[:,0]==l2),
                        np.argwhere(bLON[0,:]==l1)])
                        
    #A ZERAR AS VELOCIDADE ABAIXO DA BATIMETRIA
    RADMED = np.tile(radmed,(PVEL.shape[0],1))
    
    #APLICANDO CONDIÇÃO DE CONTORNO GVEL isobarico
    #GVEL_bar[(PVEL<RADMED)|(np.isnan(GVEL_bar))] = np.nan
    #APLICANDO CONDIÇÃO DE CONTORNO GVEL isopicnal
    #GVEL_iso[(PVEL<RADMED)|(np.isnan(GVEL_iso))] = np.nan
    #APLICANDO CONDIÇÃO DE CONTORNO VMDR
    VMDR[(PVEL<RADMED)|(np.isnan(VMDR))]         = np.nan   
    
    UMDR = -VMDR*np.sin(angle)
    VMDR = VMDR*np.cos(angle)
        
    #ENCONTRANDO PROFUNDIDADE PARA CADA PONTO DA RADIAL
    lonmed = lonctd[:-1]+np.diff(lonctd)/2
    latmed = latctd[:-1]+np.diff(latctd)/2 
    
    radmed = np.array([])
    for lonrad,latrad in zip(lonmed,latmed):
        l1 = near(bLON[0,:],lonrad,1)
        l2 = near(bLAT[:,0],latrad,1)
        radmed = np.append(radmed,BAT[np.argwhere(bLAT[:,0]==l2),
                        np.argwhere(bLON[0,:]==l1)])
    
    try:
        f = scint.interp1d(DVEL[0,:],radmed,kind='cubic')
        xbat = np.linspace(DVEL.min(),DVEL.max(),100)
        bat = f(xbat)
    except:
        f = scint.interp1d(DVEL[0,:],radmed,kind='linear')
        xbat = np.linspace(DVEL.min(),DVEL.max(),100)
        bat = f(xbat)
    
    # suavizacao inicio e fim    
    batmax = bat[-1] + np.diff(bat[-2::])
    xbat = np.hstack([-10,xbat,distmax])
    bat =  np.hstack([bat[0]+50,bat,batmax])
    
    return xbat,bat
    
#######################################################                                    
                                                                                                            
def near(dat,val,how_many=1):
    dif = np.abs(dat-val)
    idx = np.argsort(dif)
    return dat[idx][:how_many]