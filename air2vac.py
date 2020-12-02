"""
Nick Milson
Air wavelengths to vacuum wavelengths
"""

def air2vac(lambda_air):
    """
    The conversion used is from Equation (1) of Edlen (1966), Metrologia, 2, 71
    also see the review of conversions in Murphy et al. (2001), MNRAS, 327, 1223.
    Sigma is now iterated once to get better estimate of 10^4/lambda_vac.
    
    Parameters
    ----------
    lambda_air : array of floats
        input vector of air wavelengths (in A).

    Returns
    -------
    lambda_vac : array of floats
        output vector of vacuum wavelengths (A).
    """
    
    sig= 1.0e4/(lambda_air)
    term0= 8.34213e-5
    term1= 2.406030e-2/(130 - sig**2)
    term2= 1.5997e-4/(38.9 - sig**2)
    flambda= term0 + term1 + term2
    lambda_vac= lambda_air*( 1.0 + flambda )
    
    sig= 1.0e4/lambda_vac                         
    term0= 8.34213e-5
    term1= 2.406030e-2/(130 - sig**2)
    term2= 1.5997e-4/(38.9 - sig**2)
    flambda= term0 + term1 + term2
    lambda_vac= lambda_air*( 1.0e0 + flambda )
    
    return lambda_vac
