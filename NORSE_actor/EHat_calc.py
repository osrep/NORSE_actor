# EHat_calc.py
# 30/01/2019
# Written by Soma Olasz
#
# This script is created to calculate the EHat physical parameter needed to do NORSE calculations. EHat is the electric
# field normalized to the critical electric field needed for runaway electrons to be possible, derived by Connor and
# Hastie (NF, 15 (1975), 415). It will take the necessary plasma parameters, and return the value of EHat.

# import modules
import numpy as np


# Define a function to calculate the EHat variable associated with NORSE. Takes the electron number density ([m^-3]),
# the Coulomb logarithm ([na], Wesson, 'Tokamaks') and the electric field ([V/m]) as an input.
def calculate(electronDensity, CoulombLogarithm, electricField):

    # Give physical parameters
    me = 9.1e-31            # electron mass, [kg]
    c = 3e8                 # speed of light, [m/s]
    e = 1.6e-19             # elementary charge, [C]
    epsilon0 = 8.85e-12     # vacuum permittivity, [F/m]

    # Calculate critical electric field
    Ec = (electronDensity*np.power(e, 3)*CoulombLogarithm)/(4*np.pi*np.power(epsilon0, 2)*me*np.power(c, 2))    # [V/m]

    # Calculate EHat
    x = electricField/Ec

    return x
