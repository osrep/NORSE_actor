# readIN.py
# 25/01/2019
# Written by Soma Olasz
#
# This file is written to read in CPOs from the EU-IM database and give them to the NORSE script. It will take CPO
# inputs through a ualpython actor and create usable variables from them. This script only works in a Kepler workflow.


# Create a variable from the given CPO in numpy ndarray form.
# A specific element can be taken from this array using the (optional) index variable
def convert(CPO, index=None):

    # Convert the CPO into a python variable
    x = CPO[index]

    return x