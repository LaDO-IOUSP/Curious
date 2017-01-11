# -*- coding: UTF-8 -*-
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge


def pizzaplot(center, radius, angle=0,nb=2, ax=None, colors=[],**kwargs):
    ''' Plots circle with inputed number of divisions with different colors(multicolor scatter).
    ===========================================================================
    
    Input :
    
    -> center: center(x,y) of the scatter
    
    -> radius: radius of the circle (float)

    -> angle: angle of rotation of the color division (degrees [0...360])    

    -> nb: number of colors in the same plot (float)
    
    -> colors: colors to fill the scatter (list)
        
    ===========================================================================
    
    Output:
    
    -> Returns Matplotlib Patch & plots the scatter
    
    ===========================================================================
    
    Python version by:
    Hélio Almeida (helio.almeida@usp.br) 
    Dante Campagnoli Napolitano (dante.napolitano@usp.br) 
    @ LaDO-IOUSP in 11/01/2017
    
    '''
    w= []
    if len(colors)!=nb:
        raise ValueError('Number of colors and parts of scatter must be the same')
    if ax is None:
        ax = plt.gca()
    for i in np.arange(1,nb+1):
        exec('theta%s = angle+%i/%f*360.'%(str(i),i,float(nb)))
    for i in np.arange(1,nb+1):
        if i==nb:
            exec('w%s = Wedge(center, radius, theta%i, theta%i, fc=colors[%i], **kwargs)'%(str(i),i,1,i-1))
        else:
            exec('w%s = Wedge(center, radius, theta%i, theta%i, fc=colors[%i], **kwargs)'%(str(i),i,i+1,i-1))
        exec('ax.add_artist(w%i)'%i)
        exec('w.append(w%i)'%i)
    return w
