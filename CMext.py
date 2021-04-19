

"""
Interstellar extinction
"""

import numpy as np
import matplotlib.pyplot as plt

def ccmext(w):
    """
    Python function to calculate interstellar extinction curve, using 
    the formula of CCM (ApJL 329, L33, 1988). This function has been
    extended longward to 9000 A in the near IR using the revised formula
    of CCM2 (ApJ 345, 245, 1989).
    
    arguments:
    w - array of wavelengths (in A)
    
    returns:
    ax - array of a(w) coefficients to CCM2 fit
    bx = array of b(w) coefficients to CCM2 fit
    """
    
    lenw = len(w)
    ax= np.empty(lenw)
    bx= np.empty(lenw)
    xinv= (1.0e4)/w     
    for i in range(lenw):
        x = xinv[i]
        if x > 8.0: # FUV
            xx = x - 8.0
            ax[i] = -1.073 + xx*(-0.628 + xx*(0.137 - 0.070*xx))
            bx[i]= 13.670 + xx*(4.257 + xx*(-0.420 + 0.374*xx))
        elif x > 3.3:  # UV
            fa= 0.0
            fb= 0.0
            if (x >= 5.9): #correction term for FUV
                xx= x-5.9
                fa= (-0.04473 - 0.009779*xx)*(xx**2)
                fb= (0.21300 + 0.120700*xx)*(xx**2)
            ax[i]= 1.752 - 0.316*x - 0.104/((x-4.67)**2 +0.341) + fa  
            bx[i]= -3.090 + 1.825*x +1.206/((x-4.62)**2 +0.263) + fb 
        elif x > 1.1 : # optical & NIR
            y = x-1.82
            ax[i]= 1.0+y*(0.17699+y*(-0.50447+y*(-0.02427+y*(0.72085 + \
                    y*(0.01979+y*(-0.77530+y*0.32999))))))
            bx[i]= y*(1.41338+y*(2.28305+y*(1.07233+y*(-5.38434 + \
                    y*(-0.62251+y*(5.30260-y*2.09002))))))
        else: # IR - CCM expression extrapolated to infinite wavelength
            ax[i ]=  0.574*(x**1.61)
            bx[i] = -0.527*x(**1.61)
    return  ax, bx

