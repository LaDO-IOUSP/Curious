# -*- coding: utf-8 -*-
import numpy as np

def rot_vel(U,V,ang,prop='uv'):
    ''' Rotate velocities along/cross transect
===========================================================================

Input :
-> U: W-E velocity
-> V: S-N velocity
-> ang: transect angle [-180,..,180] 
-> prop: return rotated velocity ['u':along; 'v':cross; 'uv':both]
         Default prop to return is 'uv'
===========================================================================

Output:
-> Rotated velocity vectors based on input 'prop'
===========================================================================

Version 1.0.0
Hélio Almeida & Dante Napolitano @ Laboratório de Dinâmica Oceânica
Oceanographic Institute - University of São Paulo
23-09-2016
'''

    C=np.cos(np.deg2rad(ang))
    S=np.sin(np.deg2rad(ang))
    Un=U*C+V*S
    Vn=U*C-V*S
    if prop=='u':
        return Un
    elif prop=='v':
        return Vn
    elif prop=='uv':
        return Un,Vn