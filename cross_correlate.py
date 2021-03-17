def cross_correlate( wavelength, reference, spectrum, intervals, k_range = list(range(-100, 101)), 
                    M = 200, w_shift = False, graphs = False ) :
    
    """
    This is a general cross-correlation function. The idea is that the function finds the integer shift which 
    results in the highest correlation coefficient, then it fits the region around this maximum with a 
    parabola, and takes the vertex of the parabola as the best shift estimate. 
    
    By James Munday and Nick Milson


    Parameters:
        
    wavelength - an array containing the full wavelength range.

    reference - the reference spectrum to cross-correlate against.
    
    spectrum - the other spectrum being compared against reference.
    
    intervals - an array containing the start and endpoints of the intervals over which to cross-correlate. 
                For instance, if you want to cross-correlate over the intervals 6200-6400A and 6700-6800A, you 
                would pass the following: [6200,6400, 6700,6800] (obviously the spacing is for appearance).
                
    k_range - a list of the grid shift values you wish to test. The default is [-100: 101], which will check 
              shifts ranging from -100 up to 100. If you change the value to the string "variable", the 
              function will start from 0 and then move to check +/-1 and then +/-2 and continue. When it finds
              a maximum correlation coefficient, it saves the associated k-value. The function will then 
              continue searching through another M further k values (in each direction, so if you find a 
              maximum at 63 and M=100, the function will continue to check through from 63 -> 163 and from -63 
              -> -163). If it finds another maximum in this region, it resets and goes on for another M; if it
              doesn't then it takes the k value of this maximum as the best shift. To be used when the shifts 
              are very large and hard to predict.
              
    M - defines how much farther after a maximum you wish to check for another maximum. Only used when 
        k_range == "variable". Default value is 200. Large values have the advantage of checking more shifts, 
        but will slow down the program.
        
    w_shift - a Boolean variable to decide if you want the actual wavelength shift returned from the function. 
              If True, it will be returned as the third item. Default value is false.
              
    graphs - a Boolean variable for if the user would like plots made of the quadratic fit to the cross
             correlations, for every wavelength interval. False by default.
        
             
    Returns:
        
    v_shift - mean radial velocity shift
    
    cc_maxes - the cross correlation values, 
    
    d_lambda (optional) - the wavelength shift, 
             
    """
    
    # Lists for shifts, central wavelengths, and max correlation coefficient
    shift_list = []
    lambda_0 = []
    cc_maxes = []
    
    # Finding wavelength increment
    increment = wavelength[1] - wavelength[0]
    
    # Making wavelength into list
    wavelength = list( wavelength )
    
    # List for holding indices
    indices = []
    
    num_intervals = len(intervals) // 2

    # Finding the indices for the associated wavelength intervals
    for i in range( num_intervals ) :
    
        indices.append( wavelength.index( intervals[2*i] ) )
        indices.append( wavelength.index( intervals[2*i+1] ) )
    
    #--------------------------------------------------------------------------------------------------------
    # Checking if the k range should be constant or variable
    if type(k_range) == list :
        
        # Running through each interval
        for n in range( num_intervals ) :
    
            # List for correlation coefficients
            cc_list = []
    
            # Variable for the reference spectrum in the specific interval
            ref_data = np.array( reference[ indices[2*n]: indices[2*n+1] ] )
            
            # Variables for biggest correlation coefficient
            biggest = -1
            kmax = 0
            
            # List of shifts to check
            k_list = list( k_range )
            
            # Going through shifts
            for k in k_list :
                
                # Calculating correlation coefficient
                corr_coef = np.corrcoef( ref_data, spectrum[ indices[2*n]+k: indices[2*n+1]+k ] )[1, 0]
                cc_list.append(corr_coef)
                
                # Keeping track of the largest
                if corr_coef > biggest :
                    
                    biggest = corr_coef
                    kmax = k
            
            # Making sure indices don't go past edges
            if kmax in k_list[:5] :
                kmax = k_list[0] + 3
            elif kmax in k_list[-5:] :
                kmax = k_list[-1] - 3
            
            # Finding index of the maximum k value
            k_index = k_list.index(kmax)
    

            k_arr = np.linspace(k_list[0],k_list[-1],1000)

            # Getting a better estimate for shift by fitting with parabola (ax^2+bx+c)
            parameters = np.polyfit( k_list[k_index-3: k_index+4], cc_list[k_index-3: k_index+4], 2 )
            a, b,c = parameters[0], parameters[1], parameters[2]
            shift_best = -b / (2*a)
            
            if graphs:
                plt.figure()
                plt.plot(np.array(k_list)*increment,cc_list,marker = ".", ls = "")
                plt.plot(np.array(k_list[k_index-3: k_index+4])*increment,np.array(cc_list[k_index-3: k_index+4]),marker = "X", ls = "", label = "Fitted points")
                plt.plot(k_arr*increment, a*(k_arr)**2 + b * k_arr + c, label = "quadratic fit")
            
                
                plt.plot(np.ones(100)*shift_best*increment, np.linspace(min(cc_list),max(cc_list) + 0.02,100), ls = "--", label = "parabola maximum")
                plt.ylim(bottom = min(cc_list) - 0.02)
                plt.xlabel("Lag (A)")
                plt.ylabel("Cross correlation")
                plt.title(str(intervals[2*n]) + " A - "+ str(intervals[2*n+1]) + " A")
                plt.legend()
                plt.show()
            # Adding best estimate, and central wavelength onto lists
            shift_list.append( shift_best )
            
            central_wavelength = ( wavelength[ indices[2*n] ] + wavelength[ indices[2*n+1]-1 ] ) / 2
            lambda_0.append( central_wavelength  )
            
            # Adding max coefficients to list
            cc_maxes.append( max(cc_list) )
    
    
    
    #--------------------------------------------------------------------------------------------------------
    # Checking if the k range should be constant or variable
    elif k_range == "variable" :
        
        # Running through each interval
        for n in range( num_intervals ) :
            
            # List for correlation coefficients
            cc_list = []
            
            # Variable for the reference spectrum in the specific interval
            ref_data = np.array( reference[ indices[2*n]: indices[2*n+1] ] )
            
            # Variables for biggest correlation coefficient
            biggest = -1
            kmax = 0
            
            # Starting value for shift
            k = 0
            done = False
            counter = 0
            
            # Going through shifts
            while not done :
                
                # Calculating correlation coefficient for +k
                corr_coef = np.corrcoef( ref_data, spectrum[ indices[2*n]+k: indices[2*n+1]+k ] )[1, 0]
                cc_list.append(corr_coef)
                
                # Keeping track of the largest
                if corr_coef > biggest :
                    
                    biggest = corr_coef
                    kmax = k
                    counter = 0
                
                # Coefficient for -k
                if k != 0 :
                    corr_coef = np.corrcoef( ref_data, spectrum[ indices[2*n]-k: indices[2*n+1]-k ] )[1, 0]
                    cc_list = list( [corr_coef] + cc_list )
                    
                # Keeping track of the largest
                if corr_coef > biggest :
                    
                    biggest = corr_coef
                    kmax = -k
                    counter = 0
                    
                # Moving to next k
                k += 1  
                counter += 1
                
                # Checking if we've hit true maximum
                if counter > M :
                    
                    done = True
                
            # Finding index and value of the maximum k value
            k_best = kmax
            
            # Adding best estimate, and central wavelength onto lists
            shift_list.append( k_best )
            
            central_wavelength = ( wavelength[ indices[2*n] ] + wavelength[ indices[2*n+1]-1 ] ) / 2
            lambda_0.append( central_wavelength  )
            
            # Adding max coefficients to list
            cc_maxes.append( max(cc_list) )
      
    
    #===========================================================================================================
    # Converting shift values to wavelength
    d_lambda = np.array(shift_list) * increment
    
    # Calculating actual velocity shift  
    c = 2.99792458*(1e5)
    v_shift = c * ( d_lambda / lambda_0 )
    
    # Returning the actual wavelength shift if asked for it, mostly used for plotting
    if w_shift :
        return v_shift, cc_maxes, d_lambda
    
    # Returning mean velocity shift and associated error
    return v_shift, cc_maxes
