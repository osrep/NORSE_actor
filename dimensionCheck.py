# dimensionCheck.py
# 21/12/2018
# Written by Soma Olasz
#
# This function is created to check the dimensions of given input data. It will check if the given data arrays have the
# same dimensions.


# Define a function to compare the dimensions of any number of input arrays.
# This function will check if any two of the input arguments have different dimensions, and if so will raise an error.
def isIdentical(*args):

    # Go through all the input arguments.
    for i in range(0, len(args)):

        # Set a variable to identify the other array being checked.
        j = i

        # Check all the previous arrays starting from the current, ith array.
        while j != 0:

            # Lower the value of j to go though all the previous arrays.
            j = j-1

            # Check if the ith and the current jth array have the same dimensions.
            if args[i].shape != args[j].shape:

                # If not, give some information, and raise an error
                print('Input arguments', i+1, 'and', j+1, ' do not have identical dimensions.')
                print('Please check your input data!')
                raise Exception('Dimensions of two of the input quantities are not identical.')

    return
