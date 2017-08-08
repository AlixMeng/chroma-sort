import glob
import numpy as np
import matplotlib.pyplot as plt

def get_easysampler_datetime(file_name):
    """Load the Easysampler sample times and return a sorted vector of timepoints in datetime64 format"""

#    from datetime import datetime
    from pandas import to_datetime

    sample_metadata = np.genfromtxt(file_name, dtype=str,delimiter=', ')
      
    date_time = np.array(to_datetime(sample_metadata[:,1]))

    date_time = date_time[date_time.argsort()]

    return(date_time)

def get_hplc_data():
    """Get the HPLC data from .CSV files (note caps) and sort by retention time"""

    tmpfile = 'temp'
    data = []

    # Convert utf-16 to utf-8 so that np.genfromtxt can read it
    # For some reason pd.csv_read doesn't work, I would prefer it. To be fixed later.

    for filename in glob.glob('*.CSV'):
        with open(filename,'r',encoding='utf-16') as infile, open(tmpfile,'w',encoding='utf-8') as outfile:
            for line in infile:
                outfile.write(line)
    
        raw = np.genfromtxt(tmpfile,delimiter=',')

        data.append(np.array([raw[:,1],raw[:,4]]))

    # Count the rows and columns

    rows = len(data)
    columns = 0

    for i in range(len(data)):
        columns = columns + len(data[i][0])

    # Initialize the table array

    table = np.zeros((rows + 1,columns))

    # Determine the indices for inserting the individual tables into the large data table
    # Basically determining the number of columns in each individual table

    cuts=[0]

    for i in range(len(data)):
        cuts.append(len(data[i][0]))
        cuts_sum = np.cumsum(cuts)

    # Insert the individual tables into the large data table based on the indices

    for i in range(len(data)):
        table[0][cuts_sum[i]:cuts_sum[i + 1]] = data[i][0]
        table[1 + i][cuts_sum[i]:cuts_sum[i + 1]] = data[i][1]

    # Sort the table based on the retention times in row 0 with argsort()

    table_sorted = table[:,table[0].argsort()]

    return(table_sorted)

def hplc_split(data,thresh=0.3):
    """Split a sorted HPLC data array into a list of subarrays based on a retention time difference threshold"""

    # Calculate forward and backward differences. At the moment it seems that only the backward difference
    # is needed for the array splitting but maybe I will need the forward difference later.

    n = np.arange(len(data[0]))
    
    diffs_forw = np.zeros(len(data[0]))
    for i in n[:-1]:
        diffs_forw[i] = data[0][i+1] - data[0][i]

    diffs_back = np.zeros(len(data[0]))
    for i in n[1:]:
        diffs_back[i] = data[0][i] - data[0][i-1]

    # Fill the right edge of the forward difference array with a backward difference, to make an array
    # of the same length as the data array

    diffs_forw[-1] = diffs_back[-1]
    
    # Fill the left edge of the backward difference array with a forward difference
    
    diffs_back[0] = diffs_forw[0]

    # Find the indices where the difference is greater than the threshold

    ind_forw = np.ravel(np.where(diffs_forw >= thresh))
    ind_back = np.ravel(np.where(diffs_back >= thresh))

    # Create the union of the set. This was used earlier when trying to find the flat points.
    # For the moment it is not needed.

    ind_union = np.union1d(ind_forw,ind_back)

    # Split the array at the indices where the difference is greater than the threshold

    data_split = np.split(data,ind_back,axis=1)
    # Note that this cannot be converted to an array, because the rows have unequal lengths

    return(data_split)

def window_average(series, n=3):
    """Calculate the average of a (time) series, with a given window size"""
    f = len(series) - (len(series) % n) 
    # Make sure the input array to np.split is evenly divisible by window size = f/n
    
    sub = np.split(series[:f],f/n)
    
    means = np.mean(sub,axis=1)
    
    return(means)

def hplc_condense(data,thresh=3,hist=False):
    """Condense the sorted, split HPLC data, with threshold based on peak frequency"""

    # Calculate the size of each sub-array in the data array (actually a list of arrays)
    
    size = np.zeros(len(data))
    for i,sub in zip(range(len(data)),data):
        size[i] = len(sub[0])

    # Filter the sub-arrays based on the threshold
    
    peaks = np.where(size >= thresh)[0]

    # Create a new table, with the header row the average retention times, and
    # the rest of the rows as the sum over the columns -- because each row only
    # contains one measured value and the rest as zeros, this will give the
    # unique measured value in the row

    condensed = np.zeros([len(data[0]),len(peaks)],dtype=float)

    for i in range(len(peaks)):
        condensed[0,i] = np.mean(data[peaks[i]][0],axis=0)
        condensed[1:,i] = np.sum(data[peaks[i]][1:],axis=1)

    # Optionally show the sizes of each array, as a sort of histogram
        
    if(hist):  
        plt.gcf().clear()
        plt.bar(np.arange(len(data)),size,align='center')
        plt.xticks(np.arange(len(data)))
        plt.show()
        
    return(condensed)