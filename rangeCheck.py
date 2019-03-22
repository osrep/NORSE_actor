# rangeCheck.py
# 21/12/2018
# Written by Soma Olasz
#
# This function is created to check the range of given input data. It will check if the maximum or the minimum is
# within a given limit, and if not raise an exception.


# Define a function to check if the maximum or the minimum of a given data array is bigger or smaller than a given limit
# 'data' is an array
# 'minOrMax' is either 'min' or 'max'
# 'limit' is the the value of the limit to be checked
def within(data, minOrMax, limit):

    if minOrMax == 'min':

        # Check if the minimum of the given data array is smaller than the limit given
        if min(data) < limit:

            raise Exception('The given data values are not in the valid range')

    elif minOrMax == 'max':

        # Check if the maximum of the given data array is larger than the limit given
        if max(data) > limit:
            raise Exception('The given data values are not in the valid range')

    else:

        raise Exception("Please type 'min' or 'max' as the minOrMax argument!")

    return
