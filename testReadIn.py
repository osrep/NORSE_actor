# testReadIn.py
# 20/12/2018
# Written by Soma Olasz
#
# This function is created to read in the external distribution and grid vectors used in NORSE. These data are in .mat
# format at the moment. This is a TEST file and will be replaced with a function to read in CPOs in the future.

# Import modules needed for loading .mat files
import scipy.io


# Define a function for loading arbitrary variables
# 'filename' is the name of the .mat file, with extension, in string format
# 'location' is the folder of the .mat file on the computer with \ at the end, in string format
# 'variable' is the name of the variable if the file is loaded into MATLAB workspace, in string format
def load(filename, location, variable):

    # Load the variable from the given location
    matlabVariable = scipy.io.loadmat(location + filename)

    x = matlabVariable[variable]

    return x
