# pGridMode
# 02/01/2019
# Written by Soma Olasz
#
# This file is written to extract the the pGridMode of a NORSE object from a given p grid vector. This script follows
# the method written in the exampleNORSERun.m.

# Import modules needed for the function below
import numpy as np


# Define a function to extract the pGridMode
def extract(gridVector, gridMax):

    # If the difference between to adjacent coordinates is constant the grid is uniform
    if gridVector[3] - gridVector[2] == gridVector[gridMax - 1] - gridVector[gridMax - 2]:
        pGridMode = 0
        GridParameter = 0

    else:

        # Calculate the points of the uniform grid used to create the polynomial grids and the expected grid points.
        extP = gridVector[:gridMax - 1]
        extP = np.append(0, extP)
        extP = extP.reshape(gridMax, 1)

        s = np.linspace(0, 1, gridMax)
        s = s.reshape(gridMax, 1)

        # Try for each of the three possible polynomials if the external grid can be fitted. If so, save the appropriate
        # grid mode.
        for i in range(2, 5):

            GridParameter = (extP[1, 0]*np.power(s[2, 0], i) - extP[2, 0]*np.power(s[1, 0], i)) \
                            / (extP[2, 0] * s[1, 0] - extP[1, 0] * s[2, 0])
            scaleFactor = max(gridVector) / (np.power(s[-1, 0], i) + GridParameter * s[-1, 0])

            # Find an average GridParameter using the scaleFactor
            extP1 = np.delete(extP, 0, 0)    # TODO delete the zero in position 0 (can be included)
            s1 = np.delete(s, 0, 0)          # TODO delete the zero in position 0 (can be included)
            GridParameter = np.divide((extP1/scaleFactor - np.power(s1, i)), s1)
            GridParameter = np.mean(GridParameter)

            # Calculate an expected grid as NORSE will calculate it later for comparison
            expectedGrid = scaleFactor*(np.power(s, i) + GridParameter*s)
            expectedGrid1 = np.delete(expectedGrid, 0, 0)   # TODO delete the zero in position 0 (can be included)

            # Calculate the relative difference between the expected and the external grid
            diff = np.divide(abs(extP1 - expectedGrid1 ), extP1)

            # If the difference is small enough, we accept the pGridMode
            if max(diff) < 1e-12:
                pGridMode = i-1
                break

    # TODO pGridMode 4 has to be implemented

    return pGridMode, GridParameter
