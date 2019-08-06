# externalNORSERun.py
# 20/12/2018
# Written by Soma Olasz based on externalNORSERun.m
#
# This file is created to test the possibility of running NORSE code from Python by implementing the commands from the
# example file externalNORSERun.m in Python language. The example file can be found in the GitHub project for NORSE,
# in issue #4, dedicated for the modification of NORSE code to accept numerical distribution.

# Import necessary modules and Python scripts.
import readIn
import testReadIn
import rangeCheck
import dimensionCheck
import pGridMode
import matlabDouble
import EHat
import CoulombLogarithm
import numpy as np
import matlab.engine

import pickle
import ual
from time import asctime
import copy


#############################
# Load the external variables
#############################

# Save the location of the matlab variables needed for the NORSE tests
# Laptop loaction: 'D:\\ToDo\\Munka\\NORSE\\NORSE\\examples\\'
# Gateway location: '/pfs/work/g2solasz/git/NORSE/examples/'

# Load the matlab variables into Python in numpy array form
# This results in  three separate arrays
f = testReadIn.load('outputAdvanced.mat', '/pfs/work/g2solasz/git/NORSE/examples/', 'f')
extPBig = testReadIn.load('externalPBig.mat', '/pfs/work/g2solasz/git/NORSE/examples/', 'extPBig')
extXiBig = testReadIn.load('externalXiBig.mat', '/pfs/work/g2solasz/git/NORSE/examples/', 'extXiBig')

# Create a list from the external data variables
inputData = [f, extPBig, extXiBig]

##################################################
# Put some sanity checks on the external variables
##################################################

########################
# Check on element range

# TODO Should check the external distribution somehow. Checking for negative values does not work.

# Check if the values of extPBig are all positive
rangeCheck.within(inputData[1], 'min', 0)

# Check if the values of extXiBig are between -1 and 1
rangeCheck.within(inputData[2], 'min', -1)
rangeCheck.within(inputData[2], 'max', 1)

#####################
# Check on dimensions

# Check if any of the external data has different dimensions
dimensionCheck.isIdentical(inputData[0], inputData[1], inputData[2])

#############################################################
# Extract numerical parameters necessary to recreate the grid
#############################################################

# pMax, nP, and nXi parameters are needed to recreate the NORSE grid

# pMax is the maximum of the extPBig vector
extPMax = max(inputData[1])

# create a temporary array from the locations of the maximum values of the external pBig
temp = np.where(inputData[1] == inputData[1].max())

# nXi is the number of maximums in extPBig
nXi = len(temp[0])

# create a temporary array from the locations of the maximum values of the external xiBig
temp = np.where(inputData[2] == inputData[2].max())

# nP is the number of maximums+1 in extXiBig
nP = len(temp[0])+1

#################################################################
# Check if the dimensions of the external distribution is correct
#################################################################

if len(inputData[0]) != ((nP-1)*nXi+1):
    raise Exception('The dimensions of the external distribution are incorrect')

####################################################################
# Extract pGridMode, pGridParameter and xiGridMode from grid vectors
####################################################################

# pGridMode
pGrid = pGridMode.extract(inputData[1], nP)
pGridMode = pGrid[0]
pGridParameter = float(pGrid[1])

# xiGridMode
# TODO find a way to implement a similar method as in pGridMode

###########################
# Resolution and parameters
###########################
# From here, the code will follow the AdvancedNORSERun example file found in the Github project named in the description
# of this code. The section name is taken from there.

# Set constant physical parameters and numerical parameters
# TODO implement the calculation of EHat
# TODO the second variable int the readIn.convert operation is the radial coordinate (?)
# Constant physical parameters
T = readIn.convert(coreprof0[0].te.value, 0)         									# eV
n = readIn.convert(coreprof0[0].ne.value, 0)         									# m^{-3}
EHat = EHat.calculate(n,CoulombLogarithm.calculate(n,T),readIn.convert(coreprof0[0].profiles1d.eparallel.value, 0))	# E/E_c
Z = readIn.convert(coreprof0[0].profiles1d.zeff.value, 0)       							# Z_eff
B = readIn.convert(coreprof0[0].toroid_field.b0)                							# T

# Numerical parameters. Some has been taken from SimpleNORSERun.m to get sensible results
# nP = 175       # these have been determined earlier
# nXi = 35       #             -- || --
# yMax = 14      # Thermal momenta (gamma v/v_th); Don't need this as time dependent parameters are excluded
pMax = float(extPMax)
nL = 7
dt = 9e-5
tMax = 0.018
nSaveSteps = 30  # 0, save the distribution at all timesteps

##############
# Set up NORSE
##############

# Start a Matlab engine to be able to call Matlab scripts
eng = matlab.engine.start_matlab()

# Save the location of the  matlab scripts necessary for the run
# Laptop: 'D:\\ToDo\\Munka\\NORSE\\NORSE\\src'
# Gateway: '/pfs/work/g2solasz/git/NORSE/src'
# Gateway: '/pfs/work/g2solasz/git/NORSE_actor'

# Add the location of NORSE files to the Matlab path
eng.addpath('/pfs/work/g2solasz/git/NORSE/src')
eng.addpath('/pfs/work/g2solasz/git/NORSE_actor')

# Initialize an empty NORSE object
o = eng.NORSE()

# Change some settings (see NORSE.m for a complete list)
eng.setfield(o, 'nSaveSteps', nSaveSteps)
eng.setfield(o, 'includeHeatSink', 1)			# TODO we will use something different in ETS
eng.setfield(o, 'enforceStrictHeatConservation', 1)
eng.setfield(o, 'show1DTimeEvolution', 0)               # TODO we will not need this feature in ETS

# Setting the parameters to NORSE object
# All the variables must be Python float, so Matlab gets them as double. The calculation doesn't work with integers.
eng.SetParameters(o, float(nP), float(nXi), float(nL), float(pMax), float(dt), float(tMax), float(T), float(n),
                  float(Z), float(EHat), float(B), nargout=0)

# Set the grid parameters calculated earlier
eng.setfield(o, 'pGridMode', pGridMode)
# TODO The grid parameter is not exactly the same as the one which created the external grid. Has to set a limit on the error
# later.
eng.setfield(o, 'pGridParameter', pGridParameter)

# Before performing the calculation set the initialDistribution property
# to 4, corresponding to external distribution input
eng.setfield(o, 'initialDistribution', 4)

# Run NORSE in silent mode so no information is printed
eng.setfield(o, 'silent', True)

#####################
# Run the calculation
#####################

# Convert the numpy arrays into matlab doubles so the PerformCalculation method can use them
f1 = matlabDouble.convert(inputData[0])
extPBig1 = matlabDouble.convert(inputData[1])
extXiBig1 = matlabDouble.convert(inputData[2])

# Create a matlab structure from the input data given in Matlab doubles
input_structure = eng.createStructure(f1, 'f', extPBig1, 'extPBig', extXiBig1, 'extXiBig')

# Perform calculation
eng.PerformCalculation(o, input_structure, nargout=0)

#####################
# Writing output data
#####################

# Take the lastcolumn of o.f. This is the distribution function on the p-xi grid for the last time step.
temp = np.array(eng.extractDistribution(o))
# print(temp)

# Write data to given CPO
# TODO this is just a test yet, I would like to see, if writing works
# Give run and shot numbers
shot = 28906
run = 43

ntime = 1

itmp = ual.itm(shot, run)
itmp.create()

itmp.equilibriumArray.resize(ntime)
for i in range (ntime):
	itmp.equilibriumArray.array[i].eqgeometry.source = 'example'
	itmp.equilibriumArray.array[i].time = i

itmp.equilibriumArray.put()
itmp.close()