# hdf5Write.py
# 12/12/2019
# Written by Soma Olasz
#
# This function is created to write data to hdf5 file.

import h5py
import os
import pwd


def write(shotnumber, runnumber, *args):
	
	if not len(args)%2 == 0:
		raise Exception('Please input even number of arguments!')
	
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
	for i in range(0, len(args)/2):
			
		try:
			hf.create_dataset(args[2*(i+1)-1], data=args[2*i], maxshape = (None, args[2*i].shape[1]), chunks = args[2*i].shape)
		
		except RuntimeError:
			
			if args[2*i].shape[1] == hf[args[2*(i+1)-1]].shape[1]:
		
				hf[args[2*(i+1)-1]].resize(hf[args[2*(i+1)-1]].shape[0] + args[2*i].shape[0], axis = 0)			
				hf[args[2*(i+1)-1]][-args[2*i].shape[0]:] = args[2*i]
				
			else:
				print("Length of input data is different from the length of data in the hdf5 file. Please write to a different hdf5 file.")
				break
			
		
	hf.close()

	return
