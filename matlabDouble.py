# matlabDouble.py
# 08/01/2019
# Written by Soma Olasz
#
# This function is created to convert an arbitrary numpy array into Matlab double format.

import matlab.engine
import array as arr


def convert(npArray):

    # Convert the input numpy array into regular array.
    x = arr.array('d', npArray)

    # Convert this array into a list
    x = list(x)

    # Convert te list into a Matlab double
    x = matlab.double(x)

    return x
