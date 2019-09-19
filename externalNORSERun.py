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

import itertools
import pickle
import ual
from time import asctime
import copy


# Load the external variables

# Save the location of the matlab variables needed for the NORSE tests
# Laptop loaction: 'D:\\ToDo\\Munka\\NORSE\\NORSE\\examples\\'
# Gateway location: '/pfs/work/g2solasz/git/NORSE/examples/'

# Load the matlab variables into Python in numpy array form
f = testReadIn.load('outputAdvanced.mat', '/pfs/work/g2solasz/git/NORSE/examples/', 'f')
extPBig = testReadIn.load('externalPBig.mat', '/pfs/work/g2solasz/git/NORSE/examples/', 'extPBig')
extXiBig = testReadIn.load('externalXiBig.mat', '/pfs/work/g2solasz/git/NORSE/examples/', 'extXiBig')

# Create a list from the external data variables
inputData = [f, extPBig, extXiBig]

# Put some sanity checks on the external variables

# Check on element range

# TODO Should check the external distribution somehow. Checking for negative values does not work.

# Check if the values of extPBig are all positive
rangeCheck.within(inputData[1], 'min', 0)

# Check if the values of extXiBig are between -1 and 1
rangeCheck.within(inputData[2], 'min', -1)
rangeCheck.within(inputData[2], 'max', 1)

# Check on dimensions

# Check if any of the external data has different dimensions
dimensionCheck.isIdentical(inputData[0], inputData[1], inputData[2])

# Extract numerical parameters necessary to recreate the grid
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

# Check if the dimensions of the external distribution is correct

if len(inputData[0]) != ((nP-1)*nXi+1):
    raise Exception('The dimensions of the external distribution are incorrect')

# Extract pGridMode, pGridParameter and xiGridMode from grid vectors

# pGridMode
pGrid = pGridMode.extract(inputData[1], nP)
pGridMode = pGrid[0]
pGridParameter = float(pGrid[1])

# xiGridMode
# TODO find a way to implement a similar method as in pGridMode

# Resolution and parameters
# From here, the code will follow the AdvancedNORSERun example file found in the Github project named in the description
# of this code. The section name is taken from there.

# Set constant physical parameters and numerical parameters
# TODO the second variable int the readIn.convert operation is the radial coordinate in the CPOs. Using this, the calculation can be done through all radial points.

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

# Set up NORSE

# Start a Matlab engine to be able to call Matlab scripts
eng = matlab.engine.start_matlab()

# Save the location of the  matlab scripts necessary for the run
# Laptop: 'D:\\ToDo\\Munka\\NORSE\\NORSE\\src'
# Gateway: '/pfs/work/g2solasz/git/NORSE/src'
# Gateway: '/pfs/work/g2solasz/git/NORSE_actor'

# Add the location of NORSE files to the Matlab path
eng.addpath('/pfs/work/g2solasz/git/NORSE_hoppe/NORSE/src')
eng.addpath('/pfs/work/g2solasz/git/NORSE_actor')

# Initialize an empty NORSE object
o = eng.NORSE()

# Change some settings (see NORSE.m for a complete list)
eng.setfield(o, 'nSaveSteps', nSaveSteps)
eng.setfield(o, 'includeHeatSink', 1)			# TODO we will use something different in ETS
eng.setfield(o, 'enforceStrictHeatConservation', 1)
eng.setfield(o, 'show1DTimeEvolution', 0)
eng.setfield(o, 'conservativeParticleSource',1)

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

# Run the calculation

# Convert the numpy arrays into matlab doubles so the PerformCalculation method can use them
f1 = matlabDouble.convert(inputData[0])
extPBig1 = matlabDouble.convert(inputData[1])
extXiBig1 = matlabDouble.convert(inputData[2])

# Create a matlab structure from the input data given in Matlab doubles
input_structure = eng.createStructure(f1, 'f', extPBig1, 'extPBig', extXiBig1, 'extXiBig')

# Perform calculation
eng.PerformCalculation(o, input_structure, nargout=0)

# Writing output data to CPOs

# Take the data from the NORSE object, wchich will go into the CPO.
distribution = np.array(eng.extractDistribution(o)).tolist()
pBig = np.array(eng.extractPBig(o)).tolist()
xiBig = np.array(eng.extractXiBig(o)).tolist()

# flatten the data to python list, so it can be given to the CPO
distribution = list(itertools.chain.from_iterable(distribution))
pBig = list(itertools.chain.from_iterable(pBig))
xiBig = list(itertools.chain.from_iterable(xiBig))

# Write data to given CPO

# Give run and shot numbers
shot = parameters["shotnumber"]
run = parameters["run_out"]

ntime = 1	# TODO only calculate for one time step at the moment
nRho = 0	# TODO only calculate for one rho value at the moment. This will be the same as the range in the for loop for physical 		  parameters above.

itmp = ual.itm(shot, run)
itmp.create()

itmp.distributionArray.resize(ntime)

for i in range (ntime):
	
	# initialize the CPO
	itmp.distributionArray.array[i].distri_vec.resize(1)
	itmp.distributionArray.array[i].distri_vec[0].dist_func.f_expansion.resize(1)
	itmp.distributionArray.array[i].distri_vec[0].dist_func.f_expansion[0].grid.spaces.resize(3)
	
	# fill the coordinates
	for j in range (0,3):
		itmp.distributionArray.array[i].distri_vec[0].dist_func.f_expansion[0].grid.spaces[j].objects.resize(1)
		itmp.distributionArray.array[i].distri_vec[0].dist_func.f_expansion[0].grid.spaces[j].coordtype.resize(1,1)
		
		# p coordinate
		if j == 0: 
			itmp.distributionArray.array[i].distri_vec[0].dist_func.f_expansion[0].grid.spaces[j].objects[0].geo.resize((nP-1)*nXi+1,1,1,1)
			itmp.distributionArray.array[i].distri_vec[0].dist_func.f_expansion[0].grid.spaces[j].coordtype[0,0] = 123
			itmp.distributionArray.array[i].distri_vec[0].dist_func.f_expansion[0].grid.spaces[j].objects[0].geo[:,0,0,0] = pBig
			
		# xi coordinate
		elif j == 1:
			itmp.distributionArray.array[i].distri_vec[0].dist_func.f_expansion[0].grid.spaces[j].objects[0].geo.resize((nP-1)*nXi+1,1,1,1)
			itmp.distributionArray.array[i].distri_vec[0].dist_func.f_expansion[0].grid.spaces[j].coordtype[0,0] = 126
			itmp.distributionArray.array[i].distri_vec[0].dist_func.f_expansion[0].grid.spaces[j].objects[0].geo[:,0,0,0] = xiBig
			
		# rho coordinate
		else:
			itmp.distributionArray.array[i].distri_vec[0].dist_func.f_expansion[0].grid.spaces[j].objects[0].geo.resize(nRho,1,1,1)
			# 107 is the coordinate convention for rho_tor (see https://portal.eufus.eu/documentation/ITM/html/itm_enum_types__coordinate_identifier.html#itm_enum_types__coordinate_identifier). Might have to change this later.
			itmp.distributionArray.array[i].distri_vec[0].dist_func.f_expansion[0].grid.spaces[j].coordtype[0,0] = 107
			itmp.distributionArray.array[i].distri_vec[0].dist_func.f_expansion[0].grid.spaces[j].objects[0].geo[:,0,0,0] = 0.5
			
	# Write the distribution to the CPO
	itmp.distributionArray.array[i].distri_vec[0].dist_func.f_expansion[0].values.scalar.resize((nP-1)*nXi+1)
	itmp.distributionArray.array[i].distri_vec[0].dist_func.f_expansion[0].values.scalar[:] = distribution
		
	# Write the time
	itmp.distributionArray.array[i].time = i
		
itmp.distributionArray.put()
itmp.close()