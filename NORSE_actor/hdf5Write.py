# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 12:00:00 2019

@author: Soma Olasz
"""

# This function is created to write data to hdf5 file.

import h5py
import os
import pwd


def write_params(shotnumber, runnumber, input_list):
	
	# Get the location of the hdf5 output
	try:
		hdf5_base = os.environ['HDF5_BASE']
	
	except KeyError:
		
		try:
			hdf5_base = os.environ['HOME']
		
		except KeyError:
			
			hdf5_base = pwd.getpwuid(os.getuid()).pw_dir
	
	# Get the shot and run numbers into formatted strings
	str_shotnumber = str(shotnumber)
	str_runnumber = str(runnumber).zfill(4)

	# Define hdf5 output file
	location = hdf5_base + '/euitm_' + str_shotnumber + str_runnumber + '_NORSE.h5'
	
	# Init hdf5 file
	hf = h5py.File(location, 'a')
	
	# Go through all the input arguments.
	for element in input_list:
			
		try:
			hf.create_dataset(element['Name'], data=element['Data'], maxshape = (None, element['Data'].shape[1]), chunks = element['Data'].shape)
		
		except RuntimeError:
			
			if element['Data'].shape[1] == hf[element['Name']].shape[1]:
		
				hf[element['Name']].resize(hf[element['Name']].shape[0] + element['Data'].shape[0], axis = 0)			
				hf[element['Name']][-element['Data'].shape[0]:] = element['Data']
				
			else:
				print("Length of input data is different from the length of data in the hdf5 file. Please write to a different hdf5 file.")
				break
			
		
	hf.close()

	return

def write_dist(shotnumber, runnumber, input_list):
	
	# Get the location of the hdf5 output
	try:
		hdf5_base = os.environ['HDF5_BASE']
	
	except KeyError:
		
		try:
			hdf5_base = os.environ['HOME']
		
		except KeyError:
			
			hdf5_base = pwd.getpwuid(os.getuid()).pw_dir
	
	# Get the shot and run numbers into formatted strings
	str_shotnumber = str(shotnumber)
	str_runnumber = str(runnumber).zfill(4)

	# Define hdf5 output file
	location = hdf5_base + '/euitm_' + str_shotnumber + str_runnumber + '_NORSE.h5'
	
	# Init hdf5 file
	hf = h5py.File(location, 'a')
	
	# Go through all the input arguments.
	for element in input_list:
			
		try:
			hf.create_dataset(element['Name'], data=element['Data'], maxshape = (None, element['Data'].shape[1], element['Data'].shape[2]), chunks = element['Data'].shape)
		
		except RuntimeError:
			
			if element['Data'].shape[1] == hf[element['Name']].shape[1] and element['Data'].shape[2] == hf[element['Name']].shape[2]:
		
				hf[element['Name']].resize(hf[element['Name']].shape[0] + element['Data'].shape[0], axis = 0)			
				hf[element['Name']][-element['Data'].shape[0]:] = element['Data']
				
			else:
				print("Length of input data is different from the length of data in the hdf5 file. Please write to a different hdf5 file.")
				break
			
		
	hf.close()

	return

def write_RunawayRegion(shotnumber, runnumber, input_list):
	
	# Get the location of the hdf5 output
	try:
		hdf5_base = os.environ['HDF5_BASE']
	
	except KeyError:
		
		try:
			hdf5_base = os.environ['HOME']
		
		except KeyError:
			
			hdf5_base = pwd.getpwuid(os.getuid()).pw_dir
	
	# Get the shot and run numbers into formatted strings
	str_shotnumber = str(shotnumber)
	str_runnumber = str(runnumber).zfill(4)

	# Define hdf5 output file
	location = hdf5_base + '/euitm_' + str_shotnumber + str_runnumber + '_NORSE.h5'
	
	# Init hdf5 file
	hf = h5py.File(location, 'a')
	
	try:
		g1 = hf.create_group('runaway_region')
	
	except ValueError:
		
		g1 = hf.get('runaway_region')
	
	# Go through all the input arguments.
	for element in input_list:
			
		try:
			g1.create_dataset(element['Name'], data=element['Data'], maxshape = (None, element['Data'].shape[1]), chunks = element['Data'].shape)
			
		
		except RuntimeError:
			
			if element['Data'].shape[1] == g1[element['Name']].shape[1]:
		
				g1[element['Name']].resize(g1[element['Name']].shape[0] + element['Data'].shape[0], axis = 0)			
				g1[element['Name']][-element['Data'].shape[0]:] = element['Data']
				
			else:
				print("Length of input data is different from the length of data in the hdf5 file. Please write to a different hdf5 file.")
				break
			
		
	hf.close()

	return