	# Import the running_stats function
from helper import running_stats
# Write your median_bins_fits and median_approx_fits here:

import numpy as np
from astropy.io import fits

def median_bins_fits(list_of_files, B):
    
    ##@@@@@@@@ NB I've blanked out the actual bins (costly)
    
    # Get & numericise the data
    new_list_files = []
    for file in list_of_files:
        hdulist = fits.open(file)
        data = hdulist[0].data
        new_list_files.append(data)
        
    # Definitions for iterating over pixels
    n_rows = new_list_files[0].shape[0]
    n_cols = new_list_files[0].shape[1]
  
    # Vectorized mean, stdev and bounds
    mean_ = running_stats(list_of_files)[0] #np.mean(values)
    stdev = running_stats(list_of_files)[1] #np.std(values)
    l_bound = mean_ - stdev
    u_bound = mean_ + stdev
  
    # Vectorized "too-small" count
    ## fail (returns 1, shd be 0): too_small = np.where(values < l_bound)
    too_small = sum(new_list_files < l_bound) 
    
    # Pixel loop to calculate bin counts for each pixel

    width = stdev * 2 / B # this is still vectorized
    
    #@@@@@@@@@@@ bins_array = np.zeros((n_rows,n_cols,B), dtype=object)
    count_array = np.zeros((n_rows, n_cols,B), dtype=object)
    
    # Start pixel loop
    for row in range(n_rows): # 0,1,2, ..., n_rows-1
        for col in range(n_cols):

            counts_list = [] # the required counts (1-D list)
            temp_file_list = new_list_files # VECTOR to preserve original (costly)
            #@@@@@@@@@   bins=[] # the actual elements go here
            for i in range(1, B+1):
                inner_list=[] # the actual bin, contains elements (update array with this)
                local_count=0 # bin size
                for image in temp_file_list:
                    w = width[row,col]
                    lb = l_bound[row,col]+(i-1)*w
                    ub = lb+w
                    val = image[row,col]
                    if lb <= val < ub: # if that pixel is in interval...
                        local_count += 1 # increase bin size
                        count_array[row,col,i-1]+=1
                        #@@@@@@@@@@   inner_list.append(val) # add el to bin
                # temp_values.remove(val) this causes skipping! use copy.copy...

                counts_list.append(local_count) # this is know the answer but a list
                #@@@@@@@@@@@@@@@   bins.append(inner_list) 
            # Exit bin-size loop

                # Update original pixel map
                #@@@@@@@@@@@@@@@bins_array[row,col,i-1] = bins @@@@@@@@@@@
                #bins_array[row,col,i-1] = np.asarray(inner_list) #exceeds CPU limit ##recent
                #counts_list = np.asarray(counts_list)
                #############count_array[row,col,i-1] = counts_list[i-1] 
                #count_array[row,col,i-1] = np.asarray(counts_list) #exceeds CPU limit

        # End pixel loop

    
    # Convert both to array -  redundant following code edits above
    #count_array = np.asarray(count_array)
    #bins = np.asarray(bins)
    
    # Outputs
    return mean_, stdev, too_small, count_array

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def median_approx_fits(list_of_files, B):
    mean_, stdev, too_small, count_array = median_bins_fits(list_of_files, B)
    l_bound = mean_ - stdev
    u_bound = mean_ + stdev # array
    width = stdev * 2 / B # array
    N=len(list_of_files)
    left_count = too_small # array
    bin_count = 0
    
    # Define answer array - let's do this entirely element-wise
    medians_array = np.zeros_like(mean_)
    
    for row in range(mean_.shape[0]):
        for col in range(mean_.shape[1]):
            rolling_count = left_count[row,col]
            bin_index = 0
            for i in range(B): # recall this is 0,1,2,...,B-1
                if rolling_count >= (N+1)/2:
                    break
                rolling_count += count_array[row,col,i]
                bin_index +=1 # earliest possible is 0
            LB = l_bound[row,col] + (bin_index - 1) * width[row,col]
            UB = l_bound[row,col] + (bin_index) * width[row,col]
            median = (LB + UB) / 2
            
            medians_array[row,col] = median
    
    return medians_array
    
    