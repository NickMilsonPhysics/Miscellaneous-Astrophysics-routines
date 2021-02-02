import numpy as np

def lsq_bisect(x,y,wx,wy):
    """
    Parameters
    ----------
    x : x points
    y : y points
    wx : x weights
    wy : y weights
    
    Returns
    -------
    cov: coefficient vector of the bisector fit line
    """ 
    

   
    cof1= np.polyfit( x, y, 1, w = wy)
    cof2= np.polyfit( y, x, 1, w = wx)
     
    b1= cof1[0] #slope of normal fit:  y= a1 + b1*x
    b2= 1.0/cof2[0] #slope of inverted fit: x= a2 + b2*y --> y= -a2/b2 + x/b2
     
    # now compute the slope of the bisector of these two OLS lines
     
    t1= np.sqrt( 1.0 + b1**2 )
    t2= np.sqrt( 1.0 + b2**2 )
     
    b = ( b2*t1 + b1*t2 )/( t1 + t2 )     #slope of bisector line
      
     # the regression line always passes through the centroid (mean(x),mean(y)) and 
     # so we use this to find the constant term: a = mean(y) - bs*mean(x)
    
    xm= np.mean(x)
    ym= np.mean(y)
    a = ym - b*xm
    cof= np.array([b, a])  #define the coefficient vector of the bisector fit line

    return cof 
