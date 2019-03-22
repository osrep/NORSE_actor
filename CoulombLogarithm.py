# CoulombLogarithm.py
# 30/01/2019
# Written by Soma Olasz
#
# This script is written in order to calculate the Coulomb logarithm from certain plasma parameters. As a source,
# John Wesson, 'Tokamaks', Second Edition, Clarendon Press, 1997 was used.

# import modules
import numpy as np


# Define a function to calculate the Coulomb logarithm given the electron density ([m^-3]) and the electron temperature
# ([eV]) of the plasma.
def calculate(electronDensity, electronTemperature):

    # Calculate the Coulomb logarithm
    x = 14.9-0.5*np.log(electronDensity*1e-20)+np.log(electronTemperature*1e-3)       # [na]

    return x


print(calculate(4.625e19, 6026))
