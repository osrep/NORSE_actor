# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 10:07:00 2021

@author: Soma Olasz
"""

# Import necessary modules and Python scripts.
import rangeCheck
import dimensionCheck
import pGridMode
import matlabDouble
import EHat_calc
import CoulombLogarithm
import hdf5Write

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

# Start a Matlab engine to be able to call Matlab scripts
eng = matlab.engine.start_matlab()

# Test input distribution
try:
	fBig = distribution0.array[0].distri_vec[1].dist_func.f_expansion[0].values.scalar[:]
	extPBig = distribution0.array[0].distri_vec[1].dist_func.f_expansion[0].grid.spaces[0].objects[0].geo[:,0,0,0]
	extXiBig = distribution0.array[0].distri_vec[1].dist_func.f_expansion[0].grid.spaces[1].objects[0].geo[:,0,0,0]
	ext_rho_tor = distribution0.array[0].distri_vec[1].dist_func.f_expansion[0].grid.spaces[2].objects[0].geo[:,0,0,0]
	
	externalSwitch = 1
	
except IndexError:
	
	externalSwitch = 0

if externalSwitch == 1:

	ext_rho_size = ext_rho_tor.size

	# Extract numerical parameters necessary to recreate the grid
	# pMax, nP, and nXi parameters are needed to recreate the NORSE grid

	# pMax is the maximum of the extPBig vector
	extPMax = max(extPBig)

	# create a temporary array from the locations of the maximum values of the external pBig
	temp = np.where(extPBig == extPBig.max())

	# nXi is the number of maximums in extPBig
	nXi = len(temp[0])
	
	# create a temporary array from the locations of the maximum values of the external xiBig
	temp = np.where(extXiBig == extXiBig.max())

	# nP is the number of maximums+1 in extXiBig
	#nP = len(temp[0])+1
	nP = (len(fBig)/ext_rho_size -1)/nXi + 1 # TODO temporary solution for nP calculation

	# Define the variable (nP-1)*nXi+1
	coord_size = (nP-1)*nXi+1

	# Check if the dimensions of the external distribution is correct
	if len(fBig)/ext_rho_size != coord_size:
		raise Exception('The dimensions of the external distribution are incorrect')

	# Extract pGridMode, pGridParameter and xiGridMode from grid vectors

	# pGridMode
	pGrid = pGridMode.extract(extPBig, nP)
	pGridMode = pGrid[0]
	pGridParameter = float(pGrid[1])


	# Extract the distribution for the different rho coordinates
	f = np.zeros((ext_rho_size,coord_size))

	for i in range (0,ext_rho_size):
	
		f[i,:] = fBig[0:coord_size]
	
		fBig = fBig[coord_size:]


	# xiGridMode
	# TODO find a way to implement a similar method as in pGridMode
	
	# Add the variables to the matlab workspace
	eng.workspace["pMax"] = float(extPMax)
	eng.workspace["nP"] = float(nP)
	eng.workspace["nXi"] = float(nXi)
	eng.workspace["pGridMode"] = float(pGridMode)
	eng.workspace["pGridParameter"] = float(pGridParameter)
	
else:
	# make arbitrary grid parameters
	extPMax = 7				# set temporary value for pMax
	nP = 500
	nXi = 35
	
	# add the varaibles to the matlab worspace
	eng.workspace["nP"] = float(nP)
	eng.workspace["nXi"] = float(nXi)

nL = 15
eng.workspace["nL"] = float(nL)

t_in = parameters["dt_in"]
stop = parameters["stop"]
eng.workspace["t_in"] = float(t_in)
eng.workspace["stop"] = float(stop)


nSaveSteps = 0
eng.workspace["nSaveSteps"] = float(nSaveSteps)

# Set up NORSE

# Add the location of NORSE files to the Matlab path
eng.addpath('/pfs/work/g2solasz/git/NORSE_actor/ext/NORSE/src')
eng.addpath('/pfs/work/g2solasz/git/NORSE_actor/NORSE_actor')

# Initialize an empty NORSE object
eng.eval("o = NORSE();", nargout = 0)

# Change some settings (see NORSE.m for a complete list)
eng.eval("o.nSaveSteps = nSaveSteps;", nargout = 0)
eng.eval("o.includeHeatSink = 1;", nargout = 0)			# TODO we will use something different in ETS
eng.eval("o.enforceStrictHeatConservation = 1;", nargout = 0)
eng.eval("o.show1DTimeEvolution = 0;", nargout = 0)
eng.eval("o.conservativeParticleSource = 1;", nargout = 0)
eng.eval("o.runawayRegionMode = 2;", nargout = 0)
eng.eval("o.timeAdvanceMode = 0;", nargout = 0)

# Setting the parameters to NORSE object

if externalSwitch == 1:
	# Set the grid parameters calculated earlier
	eng.eval("o.pGridMode = pGridMode;", nargout = 0)
	# TODO The grid parameter is not exactly the same as the one which created the external grid. Has to set a limit on the error
	# later.
	eng.eval("o.pGridParameter = pGridParameter;", nargout = 0)

	# Before performing the calculation set the initialDistribution propert, depending on the existence of the external data
	eng.eval("o.initialDistribution = 4;", nargout = 0)

else:
	eng.eval("o.initialDistribution = 0;", nargout = 0)
	
# Run NORSE in silent mode so no information is printed
#eng.eval("o.silent = True;", nargout = 0)

# Numerical parameters
# Get the number of rho coordinates
rho_size = size(coreprof0[0].rho_tor_norm)

if externalSwitch == 1:
	if ext_rho_size != rho_size:
		raise Exception('There are different number of rho coordinates in the input CPOs')

# Initialize arrays for physical parameters
temperature = np.zeros(rho_size)
density = np.zeros(rho_size)
EHat = np.zeros(rho_size)
Z_eff = np.zeros(rho_size)
B0 = np.zeros(rho_size)
rhoTor_arr = np.zeros(rho_size)
E_parallel = np.zeros(rho_size)
E_critical = np.zeros(rho_size)
time = coreprof0[0].time

# Define physical constants
c = 3e8
e = 1.6e-19
m = 9.1e-31

# Constant physical parameters
for i in range(rho_size):
	
	# Fill physics arrays with values from CPOs
	temperature[i] = coreprof0[0].te.value[i]						# eV
	density[i] = coreprof0[0].ne.value[i]						# m^{-3}
	EHat[i] = EHat_calc.calculate(density[i],CoulombLogarithm.calculate(density[i],temperature[i]), coreprof0[0].profiles1d.eparallel.value[i])					# E/E_c	
	Z_eff[i] = coreprof0[0].profiles1d.zeff.value[i]					# Z_eff
	rhoTor_arr[i] = coreprof0[0].rho_tor[i]						# m
	B0[i] = coreprof0[0].toroid_field.b0						# T
	E_parallel[i] = coreprof0[0].profiles1d.eparallel.value[i]				# V/m
	E_critical[i] = E_parallel[i]/EHat[i]

# Initialize variables for storing calculation results
totalDistribution = []		# data storage for CPOs
finalPBig = []			# data storage for CPOs
finalXiBig = []			# data storage for CPOs
growth_rate = []
runaway_density = []
runaway_current = []
region_mask = []
region_pcs = []

# Initialize numpy arrays to store distribution and coordinates for hdf5 writing
Distribution = np.zeros((1,rho_size,(nP-1)*nXi+1))
PBig = np.zeros((1,rho_size,(nP-1)*nXi+1))
XiBig = np.zeros((1,rho_size,(nP-1)*nXi+1))


for i in range (rho_size):
	
	eng.workspace["T"] = float(temperature[i])
	eng.workspace["n"] = float(density[i])
	eng.workspace["Z"] = float(Z_eff[i])
	eng.workspace["EHat"] = float(EHat[i])
	eng.workspace["B"] = float(B0[i])
	
	# Renormalize time with the collision frequency
	coll_freq = E_critical[i]*e/(m*c)
	tMax = t_in * coll_freq
	dt = tMax/50
	
	eng.workspace["tMax"] = float(tMax)
	eng.workspace["dt"] = float(dt)

	
	if externalSwitch == 0:
		
		tTotal = t_in * stop * coll_freq	# calculate the total final time of the simulation	
		pMax = tTotal*EHat[i]/0.9 + 2		# set pMax in the first step to be larger than largest p
		
		eng.workspace["pMax"] = float(pMax)

	# All the variables must be Python float, so Matlab gets them as double. The calculation doesn't work with integers.
	eng.eval("o.SetParameters(nP, nXi, nL, pMax, dt, tMax, T, n, Z, EHat, B);", nargout=0)
	
	if externalSwitch == 1:
		# Convert the numpy arrays into matlab doubles so the PerformCalculation method can use them
		f1 = matlabDouble.convert(f[i,:])		
		extPBig1 = matlabDouble.convert(extPBig)
		extXiBig1 = matlabDouble.convert(extXiBig)
		
		eng.workspace["f"] = f1
		eng.workspace["pBig"] = extPBig1
		eng.workspace["xiBig"] = extXiBig1
		
		# Initialize a Grid object from the external Bigvectors
		eng.eval("g = Grid();", nargout = 0)
		eng.eval("g.nP = nP;", nargout = 0)
		eng.eval("g.nXi =nXi;", nargout = 0)
		eng.eval("g.pGridMode = pGridMode;", nargout = 0)
		eng.eval("g.xiGridMode = 1;", nargout = 0)
		eng.eval("g.pGridParameter = pGridParameter;", nargout = 0)
		eng.eval("g.pMax= pMax;", nargout = 0)
		eng.eval("g.InitializeGrid();", nargout = 0)
		
		eng.eval("f2D = g.MapBigVectorToGrid(f);", nargout = 0)
		
		# Create a matlab structure from the input data given in Matlab doubles for the given rho coordinate
		eng.eval("input_structure = createStructure(f, 'f', pBig, 'extPBig', xiBig, 'extXiBig', g, 'g', f2D, 'f2D');", nargout = 0)

		# Perform calculation
		eng.eval("o.PerformCalculation(input_structure);", nargout=0)
		
	else:
		eng.eval("o.PerformCalculation();", nargout=0)

	# Take the data from the NORSE object, which will go into the CPO.
	distribution = np.array(eng.eval("extractDistribution(o);")).tolist()
	pBig = np.array(eng.eval("extractPBig(o);")).tolist()
	xiBig = np.array(eng.eval("extractXiBig(o);")).tolist()
	growthRate = density[i]*eng.eval("extractGrowthRate(o);")*coll_freq
	runawayDensity = density[i]*eng.eval("extractFraction(o);")
	runawayCurrent = runawayDensity * e * c * np.sign(E_parallel[i])
	
	regionMask = eng.eval("extractMask(o);")
	regionPcs = eng.eval("extractPcs(o);")
	
	# flatten the data to python list, so it can be given to the CPO
	distribution = list(itertools.chain.from_iterable(distribution))
	pBig = list(itertools.chain.from_iterable(pBig))
	xiBig = list(itertools.chain.from_iterable(xiBig))
	growth_rate.append(growthRate)
	runaway_density.append(runawayDensity)
	runaway_current.append(runawayCurrent)
	region_mask.append(regionMask)
	region_pcs.append(regionPcs)

	# Save coordinates from first calculation to write to CPO
	if i == 0:
		finalPBig = pBig
		finalXiBig = xiBig
	
	# Check if the grids are the same for the later calculations as the saved grid
	else:
		if not  finalPBig == pBig:
			raise Exception('The p grid is not the same for rho index {} as for the first index.'.format(i+1))
			
		elif not finalXiBig == xiBig:
			raise Exception('The xi grid is not the same for rho index {} as for the first index.'.format(i+1))
	
	totalDistribution += distribution
	
	# Convert data lists to numpy arrays
	distribution = np.array(distribution)
	pBig = np.array(pBig)
	xiBig = np.array(xiBig)
	
	# Save distribution and coordinates to global vairables
	Distribution[0,i,:] = distribution
	PBig[0,i,:] = pBig
	XiBig[0,i,:] = xiBig

# Convert rho coordinates to list so it can be given to CPOs
rhoTor = rhoTor_arr.tolist()

# Write calculation results to CPO
# Give run and shot numbers needed for CPO writing
shot = parameters["shotnumber"]
run = parameters["run_out"]

# Initialize CPO structure
distribution0 = ual.distribution.distributionArray()

# Set numerical parameters for CPO writing
nCoord = 3	# number of different coordinates used (p, xi, rho_tor)

# initialize the CPO
distribution0.resize(1)
distribution0.array[0].distri_vec.resize(1)
distribution0.array[0].distri_vec[0].dist_func.f_expansion.resize(1)
distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces.resize(nCoord)

# fill the coordinates
for i in range (nCoord):
	distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[i].objects.resize(1)
	distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[i].coordtype.resize(1,1)

	# p coordinate
	if i == 0: 
		distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[i].objects[0].geo.resize((nP-1)*nXi+1,1,1,1)
		distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[i].coordtype[0,0] = 123
		distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[i].objects[0].geo[:,0,0,0] = finalPBig
			
	# xi coordinate
	elif i == 1:
		distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[i].objects[0].geo.resize((nP-1)*nXi+1,1,1,1)
		distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[i].coordtype[0,0] = 126
		distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[i].objects[0].geo[:,0,0,0] = finalXiBig
			
	# rho coordinate
	else:
		distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[i].objects[0].geo.resize(rho_size,1,1,1)
		# 107 is the coordinate convention for rho_tor (see https://portal.eufus.eu/documentation/ITM/html/itm_enum_types__coordinate_identifier.html#itm_enum_types__coordinate_identifier). Might have to change this later.
		distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[i].coordtype[0,0] = 107
		distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].grid.spaces[i].objects[0].geo[:,0,0,0] = rhoTor
			
# Write the distribution to the CPO
distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].values.scalar.resize(((nP-1)*nXi+1)*rho_size)
distribution0.array[0].distri_vec[0].dist_func.f_expansion[0].values.scalar[:] = totalDistribution

distribution0.array[0].distri_vec[0].profiles_1d.state.dens.resize(rho_size)
distribution0.array[0].distri_vec[0].profiles_1d.state.dens[:] = runaway_density
		
distribution0.array[0].distri_vec[0].profiles_1d.state.current.resize(rho_size)
distribution0.array[0].distri_vec[0].profiles_1d.state.current[:] = runaway_current
		
# Write the time
distribution0.array[0].time = time

# Reshape data for hdf5 writing
temperature = temperature.reshape(1,rho_size)
density = density.reshape(1,rho_size)
EHat = EHat.reshape(1,rho_size)
Z_eff = Z_eff.reshape(1,rho_size)
B0 = B0.reshape(1,rho_size)
rhoTor = rhoTor_arr.reshape(1,rho_size)
E_parallel = E_parallel.reshape(1,rho_size)
E_critical = E_critical.reshape(1,rho_size)
time = np.array(time).reshape(1,1)
growth_rate = np.array(growth_rate).reshape(1,rho_size)
runaway_density = np.array(runaway_density).reshape(1,rho_size)
runaway_current = np.array(runaway_current).reshape(1,rho_size)
region_mask = np.array(region_mask).reshape(rho_size,(nP-1)*nXi+1)
region_pcs = np.array(region_pcs).reshape(rho_size,nXi)
Distribution = np.array([np.transpose(Distribution[0,:,:])])
PBig = np.array([np.transpose(PBig[0,:,:])])
XiBig = np.array([np.transpose(XiBig[0,:,:])])

# Create dictionaries of the hdf5 input parameters
temperature = {"Name": 'temperature', "Data": temperature}
density = {"Name": 'density', "Data":  density}
EHat = {"Name": 'EHat', "Data":  EHat}
Z_eff = {"Name": 'Z_eff', "Data":  Z_eff}
B0 = {"Name": 'B0', "Data":  B0}
rhoTor = {"Name": 'rhoTor', "Data":  rhoTor}
E_parallel = {"Name": 'E_parallel', "Data":  E_parallel}
E_critical = {"Name": 'E_critical', "Data":  E_critical}
time = {"Name": 'time', "Data":  time}
growth_rate = {"Name": 'growth_rate', "Data":  growth_rate}
runaway_density = {"Name": 'runaway_density', "Data":  runaway_density}
runaway_current = {"Name": 'runaway_current', "Data":  runaway_current}
region_mask = {"Name": 'region_mask', "Data": region_mask}
region_pcs = {"Name": 'region_pcs', "Data": region_pcs}
Distribution = {"Name": 'Distribution', "Data":  Distribution}
PBig = {"Name": 'PBig', "Data":  PBig}
XiBig = {"Name": 'XiBig', "Data":  XiBig}

# Put dictionaries into a list
hdf5_param_data = [temperature, density, EHat, Z_eff, B0, rhoTor, E_parallel, E_critical, time, growth_rate, runaway_density, runaway_current]
hdf5_dist_data = [Distribution, PBig, XiBig]
hdf5_region_data = [region_mask, region_pcs]

# Write data to hdf5 file
hdf5Write.write_params(shot, run, hdf5_param_data)
hdf5Write.write_dist(shot, run, hdf5_dist_data)
hdf5Write.write_RunawayRegion(shot, run, hdf5_region_data)
