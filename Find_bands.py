import numpy as np 

#Function finding pairs of bands split by Zeeman splitting
#Does not work with any types of crossing bands
def find_pairs(Energy, n_bands, B, g_effective_theor, zero_index):   
    pairs = []
    for i in range(n_bands):
        split_bands = 999999    #large value
        index = 0
        for j in range(n_bands):
            #looking for g factor closest to theoretical at k = 0 - we should get clear Zeeman splitting in this point
            #control = abs(  abs((Energy[int((len(Energy)-1)/2)][i] - Energy[int((len(Energy)-1)/2)][j]))*1.602e-19 / (B*9.274e-24) + g_effective_theor) #len gives number of rows
            control = abs(  abs((Energy[zero_index,i] - Energy[zero_index,j]))*1.602e-19 / (B*9.274e-24) + g_effective_theor) #len gives number of rows
            #if(control < split_bands and control < 1): # <1 to avoid pairing of bands not split by Zeeman splitting 
            if control < split_bands:
                split_bands = control 
                index = j
        pairs.append(index)
        #print(Energy[zero_index][:])
    
    return pairs

# def numerate_bands(Energy, n_bands, pairs_prev, numerators_prev, pairs, numerators, is_first_calculation):  #first try - it works only for k = 0
#     numerators = np.zeros(n_bands)
#     if is_first_calculation:
#         numerators = np.linspace(0 , n_bands + 1 , 1)
#     else:
#         for i in range(len(pairs_prev)):
#             if pairs_prev[i] != pairs[i]:   #if pairs_prev is different than pairs it means that the bands have swapped 
#                 print('afasdf')

#     return numerators
def numerate_bands(pairs):
    numerators = np.zeros(len(pairs))
    band_number = 0
    for i in range(len(pairs)):
        if pairs[i] > i:
            numerators[i] = band_number
            numerators[pairs[i]] = band_number + 1
            band_number += 2
    return numerators


def crossing(Energy, zero_index, pairs_0, dk, approx_poly_deg, n_bands):
    pairs_k = np.zeros((len(Energy), n_bands))  #input has to be integer or tuple, creating array for pairs as function of k
    pairs_k[zero_index,:] = pairs_0 

    ############################################ POSITIVE K VECTORS ############################################################
    for i in range(zero_index, len(Energy)):    #loop over k range from 0 to k_max
        if i - zero_index < approx_poly_deg: #to avoid errors at the begginig of array where we dont have desired number of points
            points_prev = i - zero_index + 1
        else:
            points_prev = approx_poly_deg
        
        if i != zero_index: #expectation from interpolation versus E(k) from SC calculation 
            for j in range(n_bands):
                difference = 99999  #large value
                for k in range(n_bands):
                    control = np.abs(Approximated_next_value[j] - Energy[i][k])
                    #print(control)
                    if control < difference:
                        difference = control 
                        index = k 
                pairs_k[i,j] = pairs_k[i-1,index]     #if bands did not swap index=j, otherwise wa swap values of crossed bands in pairs_k


        Approximated_next_value = []
        for j in range(n_bands):    #Polynomial interpolation for every subband
            Values = []
            Poly = []
            for k in range(points_prev):
                Values.append(Energy[i-k][j])   #check the indexes !!!!!!
                Poly.append(poly_func(points_prev, (i - zero_index - k)*dk))
            #Equation takes the form
            #A + Bx + ..... + Nx**points_prev = Energy[i][pairs[i][j]] (those are energies in previous points)
            #A + Bx + ..... + Nx**points_prev = Energy[i-1][pairs[i-1][j]]
            #......
            #A + Bx + ..... + Nx**points_prev = Energy[i- points_prev][pairs[i-points_prev][j]]
            Coeffs = np.linalg.solve(Poly,Values)
            Approximated_next_value.append(polynomial_approximation(Coeffs, (i - zero_index + 1)*dk))

        
    ########################################## NEGATIVE K VECTORS #####################################################

    #Approximated_next_value = []
    for i in range(zero_index, -1, -1):
        #print(i)
        if zero_index - i < approx_poly_deg:
            points_prev = zero_index - i + 1
        else:
            points_prev = approx_poly_deg


        if i != zero_index:
            for j in range(n_bands):
                difference = 99999 #large value 
                for k in range(n_bands):
                    control = np.abs(Approximated_next_value[j] - Energy[i][k])
                    if control < difference:
                        difference = control
                        index = k
                    pairs_k[i,j] = pairs_k[i+1,index]

        Approximated_next_value = []
        for j in range(n_bands):
            Values = []
            Poly = []
            for k in range(points_prev):
                Values.append(Energy[i+k][j])
                Poly.append(poly_func(points_prev, (i-zero_index + k)*dk))  #for negative k values
            Coeffs = np.linalg.solve(Poly,Values)
            Approximated_next_value.append(polynomial_approximation(Coeffs, (i - zero_index - 1)*dk))    

    return pairs_k


#Anti-crossing - we may check if dispersion relation of two bands has oppostie curvatures and both must have stationary points, but hellical gap has the same properties


def poly_func(approx_poly_deg, k):
    poly_table = []
    for i in range(approx_poly_deg):
        poly_table.append(k**i)
    return poly_table

def polynomial_approximation(Coeffs, k_next):
    E_k = 0
    for i in range(len(Coeffs)):
        E_k += Coeffs[i]*k_next**i
    return E_k 