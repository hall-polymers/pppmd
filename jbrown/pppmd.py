#!/usr/bin/env python

# PPPMD: post-pocessing polymer molecular dynamics
# Suite of tools to do post processing calculations relevant to polymers 
# Version for LMH's MD class fall 2016
#
# Modified J. Brown 2016-10-20

import sys
import numpy as np
from math import floor

# read_lammpstrj: read in a lammps trajectory
#
# Input: fname, num_frames
#  fname: filename string or 'stdin' (or a value that evaluates to false) for reading from standard in 
#  num_frames: optional number of frames to read before stopping, defaults to reading in all frames
#  skip_beginning: skip this many frames at the beginning of the dump file
#  skip_between: skip this many frames between saved frames
#
# Output: r, ir, timestep, box_bounds, id2type, id2mol, mol2ids
#  r: num_frames by num_atoms+1 by 3 array of wrapped and unscaled coordinates (indexed by frame number then atom id)
#  ir: num_frames by num_atoms+1 by 3 array of image flags
#  timestep: num_frames length array of timesteps
#  box_bounds: 3D array to store boundaries of the box, indexed by frame, x/y/z, then lower/upper
#  id2type, id2mol: num_atoms+1 length arrays to map atom id to type and molecule id (if available, may be None)
#  mol2ids: num_mols+1 length list of atom id arrays corresponding to the molecules (if available, may be None)
#
# NOTE: assumes that the number of atoms in the simulation is fixed
# NOTE: also assumes that the coordinates are wrapped, so x or xs type coordinates are allowed but not xu or xsu
#
def read_lammpstrj(fname, num_frames=float('inf'), skip_beginning=0, skip_between=0):
	# helper function to read in the header and return the timestep, number of atoms, and box boundaries
	def read_header(f):
		f.readline() # ITEM: TIMESTEP
		timestep = int(float(f.readline()))

		f.readline() # ITEM: NUMBER OF ATOMS
		num_atoms = int(float(f.readline()))

		f.readline() # ITEM: BOX BOUNDS xx yy zz
		line = f.readline()
		line = line.split()
		xlo = float(line[0])
		xhi = float(line[1])
		line = f.readline()
		line = line.split()
		ylo = float(line[0])
		yhi = float(line[1])
		line = f.readline()
		line = line.split()
		zlo = float(line[0])
		zhi = float(line[1])

		return timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi


	#allow reading from standard input
	if not fname or fname == 'stdin':
		f = sys.stdin
	else:
		f = open(fname, 'r')
	
	# read in the initial header
	frame = 0
	init_timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi = read_header(f)

	# skip the beginning frames, if requested
	for skippedframe in range(skip_beginning):
		f.readline() # ITEM: ATOMS
		# loop over the atoms lines
		for atom in range(num_atoms):
			f.readline()
		init_timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi = read_header(f)

	# preallocate arrays, if possible
	if num_frames < float('inf'):
		alloc = num_frames
		inf_frames = False
	else:
		alloc = 1
		inf_frames = True
	timestep = np.zeros(alloc, np.int) # 1D array of timesteps
	box_bounds = np.zeros([alloc,3,2], np.float) # 3D array to store boundaries of the box, indexed by frame, x/y/z, then lower/upper

	timestep[frame] = init_timestep
	box_bounds[frame][0][0] = xlo
	box_bounds[frame][0][1] = xhi
	box_bounds[frame][1][0] = ylo
	box_bounds[frame][1][1] = yhi
	box_bounds[frame][2][0] = zlo
	box_bounds[frame][2][1] = zhi
	
	# NOTE: using num_atoms+1 here so that the arrays are indexed by their LAMMPS atom id
	r = np.zeros([alloc, num_atoms+1, 3], np.float) # 3D array of x, y, z coordinates, r[frame][id][coordinate]
	ir = np.zeros([alloc, num_atoms+1, 3], np.int) # 3D array of x, y, z image flags, r[frame][id][coordinate]

	id2mol = np.zeros(num_atoms+1, np.int) # array to map from atom id to molecule id, builds this from the first frame, if available
	id2type = np.zeros(num_atoms+1, np.int) # array to map from atom id to type, builds this from the first frame, if available


	# separately do the first ATOMS section so that we can initialize things, build the id2mol and id2type arrays, and so that the main loop starts with reading in the header
	line = f.readline()
	line = line.split()
	id_index = line.index("id") - 2
	if "mol" in line:
		mol_index = line.index("mol") - 2
	else:
		mol_index = None
	if "type" in line:
		type_index = line.index("type") - 2
	else:
		type_index = None

	if "x" in line:
		scaled = False
		x_index = line.index("x") - 2
		y_index = line.index("y") - 2
		z_index = line.index("z") - 2
	elif "xs" in line:
		scaled = True
		x_index = line.index("xs") - 2
		y_index = line.index("ys") - 2
		z_index = line.index("zs") - 2
	else:
		print("ERROR: x coordinate not found in lammps trajectory", file=sys.stderr)
		return

	if "ix" in line:
		ix_index = line.index("ix") - 2
		iy_index = line.index("iy") - 2
		iz_index = line.index("iz") - 2
	else:
		print("ERROR: x image flag not found in lammps trajectory", file=sys.stderr)
		return

	# loop over the atoms lines for the first frame separately, the rest of the frames will be read in below
	for atom in range(num_atoms):
		line = f.readline()
		line = line.split()

		# get the atom id
		my_id = int(line[id_index])

		# x, y, z coordinates
		r[frame][my_id][0] = float(line[x_index])
		r[frame][my_id][1] = float(line[y_index])
		r[frame][my_id][2] = float(line[z_index])

		# unscale, if necessary
		if scaled:
			r[frame][my_id][0] = r[frame][my_id][0]*(box_bounds[frame][0][1]-box_bounds[frame][0][0]) + box_bounds[frame][0][0]
			r[frame][my_id][1] = r[frame][my_id][1]*(box_bounds[frame][1][1]-box_bounds[frame][1][0]) + box_bounds[frame][1][0]
			r[frame][my_id][2] = r[frame][my_id][2]*(box_bounds[frame][2][1]-box_bounds[frame][2][0]) + box_bounds[frame][2][0]

		# x, y, z image flags
		ir[frame][my_id][0] = int(line[ix_index])
		ir[frame][my_id][1] = int(line[iy_index])
		ir[frame][my_id][2] = int(line[iz_index])

		# if available, buidl the i2mol and id2type arrays
		if mol_index is not None:
			id2mol[my_id] = int(line[mol_index])
		if type_index is not None:
			id2type[my_id] = int(line[type_index])
			
	# build the reverse of the id2mol array
	# this is a 2D array with rows of (potentially) varying length, so nest a numpy array into a python list
	if mol_index is not None:
		num_mols = id2mol.max()	
		mol2ids = [[]]
		for molid in range(1, num_mols+1):
			mol2ids.append(np.where(id2mol==molid)[0])
	else:
		num_mols = None
		mol2ids = None

	# loop over number of num_frames frames, if num_frames is infinite, will loop over all the frames in the file
	frame = 1 # this is the frame counter for frames actually read in
	frame_attempt = 0 # this is the actual frame count in the file (not counting the ones skipped in the beginning
	while frame < num_frames:

		frame_attempt += 1

		# try to read in a new header
		try:
			my_timestep, my_num_atoms, my_xlo, my_xhi, my_ylo, my_yhi, my_zlo, my_zhi = read_header(f)
		except:
			print("WARNING: hit end of file when reading in", fname, "at frame", skip_beginning + frame_attempt, file=sys.stderr)
			break

		# skip the frame if between frames to be read in and restart the loop
		if frame_attempt%(skip_between+1) > 0:
			f.readline() # ITEM: ATOMS
			# loop over the atoms lines
			for atom in range(num_atoms):
				f.readline()
			continue

		# if we don't know how many frames to read in, have to allocate more memeory for the arrays
		if inf_frames:
			timestep = np.append(timestep, 0)
			
			box_bounds = np.concatenate( ( box_bounds, np.zeros([1,3,2],np.float) ) )

			r = np.concatenate( ( r, np.zeros([1, num_atoms+1, 3], np.float) ) )
			ir = np.concatenate( ( ir, np.zeros([1, num_atoms+1, 3], np.float) ) )
		
		# update the timestep and box size arrays
		timestep[frame] = my_timestep
		box_bounds[frame][0][0] = my_xlo
		box_bounds[frame][0][1] = my_xhi
		box_bounds[frame][1][0] = my_ylo
		box_bounds[frame][1][1] = my_yhi
		box_bounds[frame][2][0] = my_zlo
		box_bounds[frame][2][1] = my_zhi

		f.readline() # ITEM: ATOMS
		# loop over the atoms lines
		for atom in range(num_atoms):
			line = f.readline()
			line = line.split()
	
			# get the atom id
			my_id = int(line[id_index])
	
			# x, y, z coordinates
			r[frame][my_id][0] = float(line[x_index])
			r[frame][my_id][1] = float(line[y_index])
			r[frame][my_id][2] = float(line[z_index])
	
			# unscale, if necessary
			if scaled:
				r[frame][my_id][0] = r[frame][my_id][0]*(box_bounds[frame][0][1]-box_bounds[frame][0][0]) + box_bounds[frame][0][0]
				r[frame][my_id][1] = r[frame][my_id][1]*(box_bounds[frame][1][1]-box_bounds[frame][1][0]) + box_bounds[frame][1][0]
				r[frame][my_id][2] = r[frame][my_id][2]*(box_bounds[frame][2][1]-box_bounds[frame][2][0]) + box_bounds[frame][2][0]

			# x, y, z image flags
			ir[frame][my_id][0] = int(line[ix_index])
			ir[frame][my_id][1] = int(line[iy_index])
			ir[frame][my_id][2] = int(line[iz_index])
	
		frame += 1

	return r, ir, timestep, box_bounds, id2type, id2mol, mol2ids

# MSD: mean squared displacement
#
# Input: r, ir, box_bounds, id2type
# computes separate MSDs for each type as given in id2type, msd_dict["type"] is a dict of these MSDs
#  r: unscaled (but wrapped) coordinates 
#  ir: image flags
#  box_bounds: boundaries of the box
#  id2type: array that maps atom id to type 
#  (format as read in from read_lammpstrj)
#
# Output: msd_dict filled with the calculated MSDs 
#  each entry is an array indexed by frame gives average MSD of all beads of the different types
#
# NOTE: does NOT do any sort of block averaging
# NOTE: assumes mass = 1 for all beads
# NOTE: does not account for changes in box size
#
def MSD(r, ir, box_bounds, id2type=[]):
	# set up some constants
	frames = len(r)
	box_size = np.array([ box_bounds[0][0][1] - box_bounds[0][0][0], box_bounds[0][1][1] - box_bounds[0][1][0], box_bounds[0][2][1] - box_bounds[0][2][0] ])

	#  allocate an array for the box center of mass which needs to be subtracted off
	box_com = np.zeros([frames,3], np.float)

	# preallocate msd vectors
	msd_dict = {}
	for type_id in set(id2type):
		msd_dict[type_id] = np.zeros(frames, np.float)

	# loop over frames
	for t in range(frames):
		# calculate the center of mass of the entire box
		for atom in range(1, len(r[0])):
			box_com[t] += r[t][atom] + ir[t][atom]*box_size
		box_com[t] = box_com[t]/(len(r[0])-1)


		# loop over atoms
		for atom in range(1, len(id2type)):
			# calculate how much the bead has moved reletive to the center of mass (note that this is a vector equation)
			diff = (r[t][atom] + ir[t][atom]*box_size - box_com[t]) - (r[0][atom] + ir[0][atom]*box_size - box_com[0])
			# the mean squared displacement is this difference dotted with itself
			msd_dict[id2type[atom]][t] += diff.dot(diff)


	# scale MSD by the number of beads of each type, to get the average MSD
	for type_id in set(id2type):
		msd_dict[type_id] = msd_dict[type_id]/sum(id2type == type_id)
	del msd_dict[0] # this is needed since id2type has a dummy entry of 0 at index 0 so that it is indexed by LAMMPS atom_id

	return msd_dict


# gofr: radial distribution function
#
# Input: r, box_bounds, bin_size
# computes g(r) between all atoms passed in the array r
#  r: unscaled (but wrapped) coordinates 
#  box_bounds: boundaries of the box
#  (format as read in from read_lammpstrj)
#  bin_size: (float) size in units of unscaled distance to bin the distances
#
# Output: (gofr, r_out)
#  gofr: array containing the radial distribution function between all read in atoms
#  r_out: array containing the corresponding distances for the gofr array
#
# NOTE: does not account for changes in box size
# NOTE: does not distinguish between atoms on the same molecule and other atoms
#
def gofr(r, box_bounds, bin_size):
	# set up some constants
	frames = len(r)
	num_atoms = len(r[0])-1
	box_size = np.array([ box_bounds[0][0][1] - box_bounds[0][0][0], box_bounds[0][1][1] - box_bounds[0][1][0], box_bounds[0][2][1] - box_bounds[0][2][0] ])
	density_tot = num_atoms/(box_size[0]*box_size[1]*box_size[2])
	num_bins = int(np.floor(0.5*min(box_size)/bin_size))

	# preallocate output vectors
	gofr = np.zeros(num_bins, np.float)
	r_out = np.linspace(0.5*bin_size, (num_bins-0.5)*bin_size, num_bins)

	# loop over frames
	for t in range(frames):
		# loop over all atoms
		for atom_id in range(1,num_atoms):
			# now find the distance from this atom to all the following atoms 
			# using numpy to vectorize the calculation for speed
			delta_r = np.abs(r[t][atom_id+1:] - r[t][atom_id]) # creates an array of vectors containing the absolute differences in the x, y, and z directions
			delta_r = np.where(delta_r > 0.5 * box_size, delta_r - box_size, delta_r) # this does an element by element comparison, wherever the absolute difference is greater than half the box size, we should have used a periodic image of one of the two atoms, so to correct this, shift that value by the appropriate box dimension. Note that since this happens on an element by element basis, it may shift the difference in the x-direction by the size of the box in the x-direction, but leave the rest of the vector alone.
			dist = np.sqrt((delta_r ** 2).sum(axis=1)) # now calculate the distances, the axis=1 option makes the sums happen over the vectors not the entire dataset, so the resulting data is an array containg the distances from atom_id
			
			# bin the distances into the gofr array
			dist = dist/bin_size
			for n in range(len(dist)):
				bin_num = int(np.floor(dist[n]))
				if bin_num < num_bins:
					gofr[bin_num] += 1.0

	# rescale based on number of atoms, bin volume, density, and number of frames
	for bin_num in range(num_bins):
		r_inner = bin_num*bin_size
		r_outer = (bin_num+1)*bin_size
		V = (4*np.pi/3)*(r_outer**3 - r_inner**3)
		gofr[bin_num] = gofr[bin_num]/((num_atoms-1)/2.0)/V/density_tot/frames

	return (gofr, r_out)


# Sofk: static structure factor
#
# Input: r, box_bounds, num_bins
# computes S(k), S(k_x), S(k_y), and S(k_z) of the atoms passed in in array r
#  r: unscaled coordinates 
#  box_bounds: boundaries of the box
#  (format as read in from read_lammpstrj)
#  num_bins: (integer) number of bins, note that the calculation grows with num_bins^3. The maximum k value will be num_bins*2*pi/box_dimension
#
# Output: (Sofk, k_out, Sofk_x, k_out_x, Sofk_y, k_out_y, Sofk_z, k_out_z) 
#  Sofk:    array containing the static structure factor between all read in atoms
#  k_out:   array containing the corresponding k values for Sofk
#  Sofk_x:  array containing the static structure factor between all read in atoms in the x direction
#  k_out_x: array containing the corresponding k values for Sofk_x
#  Sofk_y:  array containing the static structure factor between all read in atoms in the y direction
#  k_out_y: array containing the corresponding k values for Sofk_y
#  Sofk_z:  array containing the static structure factor between all read in atoms in the z direction
#  k_out_z: array containing the corresponding k values for Sofk_z
#
# NOTE: does not account for changes in box size
# NOTE: does not distinguish between atoms on the same molecule and other atoms
#
def Sofk(r, box_bounds, num_bins):
	# set up some constants
	frames = len(r)
	num_atoms = len(r[0])-1
	box_size = np.array([ box_bounds[0][0][1] - box_bounds[0][0][0], box_bounds[0][1][1] - box_bounds[0][1][0], box_bounds[0][2][1] - box_bounds[0][2][0] ])
	dk = 2*np.pi/box_size # the spacing of the k-vecors depends on the box dimensions
	min_dk = min(dk)

	# preallocate output vectors
	Sofk = np.zeros(num_bins, np.float)
	Sofk_x = np.zeros(num_bins, np.float)
	Sofk_y = np.zeros(num_bins, np.float)
	Sofk_z = np.zeros(num_bins, np.float)
	count_k = np.zeros(num_bins, np.float)
	count_k_x = np.zeros(num_bins, np.float)
	count_k_y = np.zeros(num_bins, np.float)
	count_k_z = np.zeros(num_bins, np.float)
	k_out = np.linspace(min_dk, num_bins*min_dk, num_bins)
	k_out_x = np.linspace(dk[0], num_bins*dk[0], num_bins)
	k_out_y = np.linspace(dk[1], num_bins*dk[1], num_bins)
	k_out_z = np.linspace(dk[2], num_bins*dk[2], num_bins)

	# loop over frames
	for t in range(frames):
		# loop over a grid of k-vectors inside the sphere of all indices that makes all k-vectors within the maximun range
		for ix in range(num_bins + 1):
			k_x = ix*dk[0]
			for iy in range(int(np.floor(np.sqrt(num_bins**2 - ix**2))) + 1):
				k_y = iy*dk[1]
				for iz in range(int(np.floor(np.sqrt(num_bins**2 - ix**2 - iy**2))) + 1):
					k_z = iz*dk[2]
					k = np.sqrt(k_x**2 + k_y**2 + k_z**2)
					if k == 0.: # ignore k=0
						continue

					# determine which bin the k-vector falls into 
					# this ensures that the bin at m*dk goes from (m-1/2)*dk to (m+1/2)*dk
					bin_k = int(round(k/min_dk))  - 1 
					bin_k_x = int(round(k/dk[0])) - 1
					bin_k_y = int(round(k/dk[1])) - 1
					bin_k_z = int(round(k/dk[2])) - 1

					# numpy magic to vectorize a loop over all atoms and add up the exponential terms
					# np.inner does the dot product between the k-vector and each of the vectors in the r array
					# 1.j is the imaginary unit
					Sk_temp = abs( sum( np.exp( -1.j * np.inner( [k_x,k_y,k_z], r[t][1:]) ) ) )**2 / num_atoms

					# Now bin the S(k) we just calculated into the right vectors
					if bin_k < num_bins:
						Sofk[bin_k] += Sk_temp
						count_k[bin_k] += 1
					if bin_k_x < num_bins and np.arccos(k_x/k) < 0.1618: ## count all k-vectors within 15 degrees of the x-axis as in the x-direction to improve statistics
						Sofk_x[bin_k_x] += Sk_temp
						count_k_x[bin_k_x] += 1
					elif bin_k_y < num_bins and np.arccos(k_y/k) < 0.1618: 
						Sofk_y[bin_k_y] += Sk_temp
						count_k_y[bin_k_y] += 1
					elif bin_k_z < num_bins and np.arccos(k_z/k) < 0.1618: 
						Sofk_z[bin_k_z] += Sk_temp
						count_k_z[bin_k_z] += 1
						
	# rescale by the counts in each bin
	Sofk /= count_k
	Sofk_x /= count_k_x
	Sofk_y /= count_k_y
	Sofk_z /= count_k_z

	return (Sofk, k_out, Sofk_x, k_out_x, Sofk_y, k_out_y, Sofk_z, k_out_z)

# end2end_autocorr: end to end autocorrelation function
#
# Input: r, ir, box_bounds, mol2ids
#  r: unscaled (but wrapped) coordinates
#  ir: image flags
#  box_bounds: boundaries of the box
#  mol2ids: list of atom id arrays of the molecules to be used
#  (format as read in from read_lammpstrj)
#
# Output: end to end vector autocorrelation function (1D array indexed by frame count) constructued by doting the original end to end vector to the 
#
# NOTE: all the listings in mol2ids will be used and averaged together
# NOTE: it is assumed that the end-to-end vector is the one between the lowest and highest id in each molecule (if this is not the case, you'd have to mess with mol2ids, e.g. make it only contain the ids of the two end beads)
# NOTE: scaled by the average end-to-end vector at frame 0, so that e2e_autocorr[0]=1.0
#
def end2end_autocorr(r, ir, box_bounds, mol2ids):
	frames = len(r)
	mols = len(mol2ids)

	# preallocate e2e vector arrays and autocorr array
	e2e_t = np.zeros([mols, 3], np.float)
	e2e_0 = np.zeros([mols, 3], np.float)
	e2e_autocorr = np.zeros(frames, np.float)

	# loop over time
	for t in range(frames):
		box_size = np.array([ box_bounds[t][0][1] - box_bounds[t][0][0], box_bounds[t][1][1] - box_bounds[t][1][0], box_bounds[t][2][1] - box_bounds[t][2][0] ])

		# loop over molecules
		for molid in range(1, mols):
			# assume that the ends of the chain have the maximum and minimum id numbers
			id1 = mol2ids[molid].min()
			id2 = mol2ids[molid].max()

			# calculate the end-to-end vector
			r1 = r[t][id1] + ir[t][id1]*box_size
			r2 = r[t][id2] + ir[t][id2]*box_size

			e2e_t[molid] = r2 - r1
			if t == 0:
				e2e_0[molid] = e2e_t[molid]

			# take dot products
			e2e_autocorr[t] += np.dot(e2e_0[molid], e2e_t[molid])

	# scaling
	e2e_autocorr = e2e_autocorr/(mols-1)
	e2e_autocorr = e2e_autocorr/e2e_autocorr[0]

	return e2e_autocorr


# Ion_Dynamics
# Author:  Kevin Shen, May 2019
# special_read: read in a lammps trajectory, but only extract data of certain bead types
#
# Input: fname, num_frames
#  fname: filename string or 'stdin' (or a value that evaluates to false) for reading from standard in
#  num_frames: optional number of frames to read before stopping, defaults to reading in all frames
#  skip_beginning: skip this many frames at the beginning of the dump file
#  skip_between: skip this many frames between saved frames
#
# Output: r, ir, timestep, box_bounds, id2type, id2mol, mol2ids
#  r: num_frames by num_atoms+1 by 3 array of wrapped and unscaled coordinates (indexed by frame number then atom id)
#  ir: num_frames by num_atoms+1 by 3 array of image flags
#  timestep: num_frames length array of timesteps
#  box_bounds: 3D array to store boundaries of the box, indexed by frame, x/y/z, then lower/upper
#  id2type, id2mol: num_atoms+1 length arrays to map atom id to type and molecule id (if available, may be None)
#
# NOTE: assumes that the number of atoms in the simulation is fixed
# NOTE: also assumes that the coordinates are wrapped, so x or xs type coordinates are allowed but not xu or xsu
#
def special_read(fname, types, num_frames=float('inf'), skip_beginning=0, skip_between=0):
	print("Reading configurations...")  # JB

	# helper function to read in the header and return the timestep, number of atoms, and box boundaries
	def read_header(f):
		f.readline()  # ITEM: TIMESTEP
		timestep = int(float(f.readline()))  # JB

		f.readline()  # ITEM: NUMBER OF ATOMS
		num_atoms = int(float(f.readline()))  # JB

		f.readline()  # ITEM: BOX BOUNDS xx yy zz
		line = f.readline()
		line = line.split()
		xlo = float(line[0])
		xhi = float(line[1])
		line = f.readline()
		line = line.split()
		ylo = float(line[0])
		yhi = float(line[1])
		line = f.readline()
		line = line.split()
		zlo = float(line[0])
		zhi = float(line[1])

		return timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi

	# allow reading from standard input
	if not fname or fname == 'stdin':
		f = sys.stdin
	else:
		f = open(fname, 'r')

	# read in the initial header
	frame = 0
	init_timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi = read_header(f)

	# skip the beginning frames, if requested
	for skippedframe in range(skip_beginning):
		print("Skipping " + str(skippedframe))  # JB
		f.readline()  # ITEM: ATOMS
		# loop over the atoms lines
		for atom in range(num_atoms):
			f.readline()
		init_timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi = read_header(f)

	# preallocate arrays, if possible
	if num_frames < float('inf'):
		alloc = num_frames
		inf_frames = False
	else:
		alloc = 1
		inf_frames = True
	timestep = np.zeros(alloc, np.int)  # 1D array of timesteps
	box_bounds = np.zeros([alloc, 3, 2],
						  np.float)  # 3D array to store boundaries of the box, indexed by frame, x/y/z, then lower/upper

	timestep[frame] = init_timestep
	box_bounds[frame][0][0] = xlo
	box_bounds[frame][0][1] = xhi
	box_bounds[frame][1][0] = ylo
	box_bounds[frame][1][1] = yhi
	box_bounds[frame][2][0] = zlo
	box_bounds[frame][2][1] = zhi

	# separately do the first ATOMS section so that we can initialize things, build the id2mol and id2type arrays, and so that the main loop starts with reading in the header
	line = f.readline()
	line = line.split()
	id_index = line.index("id") - 2
	if "mol" in line:
		mol_index = line.index("mol") - 2
	else:
		mol_index = None
	if "type" in line:
		type_index = line.index("type") - 2
	else:
		type_index = None

	if "x" in line:
		scaled = False
		x_index = line.index("x") - 2
		y_index = line.index("y") - 2
		z_index = line.index("z") - 2
	elif "xs" in line:
		scaled = True
		x_index = line.index("xs") - 2
		y_index = line.index("ys") - 2
		z_index = line.index("zs") - 2
	else:
		print("ERROR: x coordinate not found in lammps trajectory", file=sys.stderr)  # JB
		return

	if "ix" in line:
		ix_index = line.index("ix") - 2
		iy_index = line.index("iy") - 2
		iz_index = line.index("iz") - 2
	else:
		print("ERROR: x image flag not found in lammps trajectory", file=sys.stderr)  # JB
		return

	# NOTE: using num_atoms+1 here so that the arrays are indexed by their LAMMPS atom id
	r_temp = np.zeros([num_atoms + 1, 3], np.float)  # 3D array of x, y, z coordinates, r[frame][id][coordinate]
	ir_temp = np.zeros([num_atoms + 1, 3], np.int)  # 3D array of x, y, z image flags, r[frame][id][coordinate]

	# id2mol = np.zeros(num_atoms+1, np.int) # array to map from atom id to molecule id, builds this from the first frame, if available
	index2type = np.zeros(num_atoms + 1,
						  np.int)  # array to map from atom id to type, builds this from the first frame, if available
	index2id = np.zeros(num_atoms + 1,
						np.int)  # array to map from atom id to type, builds this from the first frame, if available

	num_types = 0
	# loop over the atoms lines for the first frame separately, the rest of the frames will be read in below
	for atom in range(num_atoms):
		line = f.readline()
		line = line.split()

		# get the atom id
		my_id = int(line[id_index])

		# build the index2type array
		my_type = int(line[type_index])
		index2type[my_id] = my_type
		if my_type in types:
			num_types += 1

		# x, y, z coordinates
		r_temp[my_id][0] = float(line[x_index])
		r_temp[my_id][1] = float(line[y_index])
		r_temp[my_id][2] = float(line[z_index])

		# unscale, if necessary
		if scaled:
			r_temp[my_id][0] = r_temp[my_id][0] * (box_bounds[frame][0][1] - box_bounds[frame][0][0]) + \
							   box_bounds[frame][0][0]
			r_temp[my_id][1] = r_temp[my_id][1] * (box_bounds[frame][1][1] - box_bounds[frame][1][0]) + \
							   box_bounds[frame][1][0]
			r_temp[my_id][2] = r_temp[my_id][2] * (box_bounds[frame][2][1] - box_bounds[frame][2][0]) + \
							   box_bounds[frame][2][0]

		# x, y, z image flags
		ir_temp[my_id][0] = int(line[ix_index])
		ir_temp[my_id][1] = int(line[iy_index])
		ir_temp[my_id][2] = int(line[iz_index])

	# NOTE: using num_types+1 here so that the arrays are indexed by their LAMMPS atom id
	num_not_types = num_atoms - num_types
	r = np.zeros([alloc, num_types + 1, 3], np.float)  # 3D array of x, y, z coordinates, r[frame][id][coordinate]
	ir = np.zeros([alloc, num_types + 1, 3], np.int)  # 3D array of x, y, z image flags, r[frame][id][coordinate]

	id2type = np.zeros(num_types + 1,
					   np.int)  # array to map from atom id to type, builds this from the first frame, if available
	id2index = np.zeros(num_types + 1,
						np.int)  # array to map from atom id to index, builds this from the first frame, if available

	# store the temporary data into real arrays
	my_id = 0
	for atom in range(num_atoms):
		index = atom + 1
		if index2type[index] in types:
			my_id += 1
			# x, y, z coordinates
			r[frame][my_id][0] = r_temp[index][0]
			r[frame][my_id][1] = r_temp[index][1]
			r[frame][my_id][2] = r_temp[index][2]

			# x, y, z image flags
			ir[frame][my_id][0] = ir_temp[index][0]
			ir[frame][my_id][1] = ir_temp[index][1]
			ir[frame][my_id][2] = ir_temp[index][2]

			id2type[my_id] = index2type[index]
			id2index[my_id] = index
			index2id[index] = my_id

	# loop over number of num_frames frames, if num_frames is infinite, will loop over all the frames in the file
	frame = 1  # this is the frame counter for frames actually read in
	frame_attempt = 0  # this is the actual frame count in the file (not counting the ones skipped in the beginning
	while frame < num_frames:

		frame_attempt += 1

		# try to read in a new header
		try:
			my_timestep, my_num_atoms, my_xlo, my_xhi, my_ylo, my_yhi, my_zlo, my_zhi = read_header(f)
		except:
			print("WARNING: hit end of file when reading in", fname, "at frame", skip_beginning + frame_attempt,
				  file=sys.stderr)  # JB
			break

		# skip the frame if between frames to be read in and restart the loop
		if frame_attempt % (skip_between + 1) > 0:
			f.readline()  # ITEM: ATOMS
			# loop over the atoms lines
			for atom in range(num_atoms):
				f.readline()
			continue

		# if we don't know how many frames to read in, have to allocate more memeory for the arrays
		if inf_frames:
			timestep = np.append(timestep, 0)

			box_bounds = np.concatenate((box_bounds, np.zeros([1, 3, 2], np.float)))

			r = np.concatenate((r, np.zeros([1, num_types + 1, 3], np.float)))
			ir = np.concatenate((ir, np.zeros([1, num_types + 1, 3], np.float)))

		# update the timestep and box size arrays
		timestep[frame] = my_timestep
		box_bounds[frame][0][0] = my_xlo
		box_bounds[frame][0][1] = my_xhi
		box_bounds[frame][1][0] = my_ylo
		box_bounds[frame][1][1] = my_yhi
		box_bounds[frame][2][0] = my_zlo
		box_bounds[frame][2][1] = my_zhi

		f.readline()  # ITEM: ATOMS
		# loop over the atoms lines
		for atom in range(num_atoms):
			line = f.readline()
			line = line.split()

			# get the atom id
			index = int(line[id_index])
			if index2type[index] in types:
				my_id = index2id[index]

				# x, y, z coordinates
				r[frame][my_id][0] = float(line[x_index])
				r[frame][my_id][1] = float(line[y_index])
				r[frame][my_id][2] = float(line[z_index])

				# unscale, if necessary
				if scaled:
					r[frame][my_id][0] = r[frame][my_id][0] * (box_bounds[frame][0][1] - box_bounds[frame][0][0]) + \
										 box_bounds[frame][0][0]
					r[frame][my_id][1] = r[frame][my_id][1] * (box_bounds[frame][1][1] - box_bounds[frame][1][0]) + \
										 box_bounds[frame][1][0]
					r[frame][my_id][2] = r[frame][my_id][2] * (box_bounds[frame][2][1] - box_bounds[frame][2][0]) + \
										 box_bounds[frame][2][0]

				# x, y, z image flags
				ir[frame][my_id][0] = int(line[ix_index])
				ir[frame][my_id][1] = int(line[iy_index])
				ir[frame][my_id][2] = int(line[iz_index])

		print("Reading frame {}".format(frame))  # JB
		frame += 1

	print('===============Summary===============')  # JB
	print('Total number of beads =', num_atoms)  # JB
	print('Total number of selected beads =', num_types)  # JB

	return r, ir, timestep, box_bounds, id2type, id2index


def buildnlist(r, bin2id, id2bin, bins, boxsize, id2type, dist_range, nearest):
	print('Building neighbor lists...')  # JB
	# loop over the local area of each bead, and shift at the periodic BCs
	id2neighbors = [[[] for n in range(len(r[0]))] for n in range(len(r))]  # JB
	for t in range(len(r)):
		# loop over atoms
		for atom in range(1, len(r[0])):
			shift = np.zeros(3, np.float)
			binloc = id2bin[t][atom]
			if nearest:
				dmin = dist_range ** 2
				dmin_id = 0
			for i in range(binloc[0] - 1, binloc[0] + 2):
				if i == -1:
					i = bins[0] - 1
					shift[0] = -boxsize[0]
				elif i == bins[0]:
					i = 0
					shift[0] = boxsize[0]
				else:
					shift[0] = 0

				for j in range(binloc[1] - 1, binloc[1] + 2):
					if j == -1:
						j = bins[1] - 1
						shift[1] = -boxsize[1]
					elif j == bins[1]:
						j = 0
						shift[1] = boxsize[1]
					else:
						shift[1] = 0

					for k in range(binloc[2] - 1, binloc[2] + 2):
						if k == -1:
							k = bins[2] - 1
							shift[2] = -boxsize[2]
						elif k == bins[2]:
							k = 0
							shift[2] = boxsize[2]
						else:
							shift[2] = 0

						# loop over the beads in this box and calculate distance
						for test_id in bin2id[t][i][j][k]:
							if (test_id in id2neighbors[t][atom]) or (test_id == atom) or (
									id2type[test_id] == id2type[atom]):
								continue
							dr = r[t][test_id] + shift - r[t][atom]
							dr2 = dr.dot(dr)
							if dr2 < dist_range ** 2:
								if nearest:
									if dr2 < dmin:
										dmin = dr2
										dmin_id = test_id
								else:
									id2neighbors[t][atom].append(test_id)
									id2neighbors[t][test_id].append(atom)

			if nearest and dmin_id != 0:
				id2neighbors[t][atom].append(dmin_id)

	return id2neighbors


# Input: r, box_bounds
#  r: unscaled (but wrapped) coordinates
#  boxsize: box dimensions (x/y/z)
# Output: bin2id, id2bin
# NOTE: does not account for changes in box size
#
def binbox(r, boxsize, bound_lo, dist_range):  # JB
	print('Binning box...')  # JB
	bins = np.floor(boxsize / dist_range).astype(int)
	id2bin = np.zeros([len(r), len(r[0]), 3], np.int)
	bin2id = [[[[[] for k in range(bins[2])] for j in range(bins[1])] for i in range(bins[0])] for t in range(len(r))]

	# loop over frames
	for t in range(len(r)):
		# loop over atoms
		for atom in range(1, len(r[0])):
			xbin, ybin, zbin = np.floor((r[t][atom] - bound_lo) / boxsize * bins).astype(int)
			try:
				id2bin[t][atom] = (xbin, ybin, zbin)
				bin2id[t][xbin][ybin][zbin].append(atom)
			except:
				print('Warning: atom {} at timestep {} slightly outside of box: '.format(atom, t),
					  file=sys.stderr)  # JB
				print('{} > {}'.format(r[t][atom], bound_hi), file=sys.stderr)  # JB

	return bin2id, id2bin, bins


# NOTE: Based on https://lammps.sandia.gov/threads/msg32219.html, atoms can be slightly outside the periodic boundary
# if such thing happens, we wrap it back
def wrap(r, box_bounds):
	print('Wrapping...')  # JB
	bound_lo = np.array([box_bounds[0][0][0], box_bounds[0][1][0], box_bounds[0][2][0]])
	bound_hi = np.array([box_bounds[0][0][1], box_bounds[0][1][1], box_bounds[0][2][1]])
	boxsize = bound_hi - bound_lo
	for t in range(len(r)):
		for atom in range(1, len(r[0])):
			shift = np.zeros(3, np.float)
			for axis, coord in enumerate(r[t][atom]):
				if coord >= bound_hi[axis]:
					shift[axis] = -boxsize[axis]
				elif coord < bound_lo[axis]:
					shift[axis] = boxsize[axis]
			if np.sum(np.abs(shift)) != 0:
				r[t][atom] += shift
	return r, boxsize, bound_lo


def ipcorr(r, id2neighbors):
	# print 'Calculating ion pair correlation function...'
	h = np.zeros([len(r), len(r[0])], np.float)
	corr = np.zeros([len(r), len(r[0])], np.float)
	for t in range(len(r)):
		for atom in range(1, len(r[0])):
			# check if the list is empty (see if there is ion pair)
			if id2neighbors[t][atom]:
				h[t][atom] = 1
				if t == 0:
					corr[t][atom] = 1
				elif id2neighbors[t][atom] == id2neighbors[0][atom] and corr[t - 1][atom] != 0:
					corr[t][atom] = 1
	return corr.mean(axis=1) / h.mean()


def cluster_autocorr(frame, sflag):
	global frame_0, sigma_start
	sigma_current = 0.0  # numerator of ACF
	index = 0

	# establish where ions are at the beginning
	if frame == 0:
		frame_0 = []
		sigma_start = 0.0  # denominator of ACF
		for ii in range(len(sflag)):
			for jj in range(ii + 1, len(sflag)):
				if sflag[ii] == sflag[jj]:
					sigma_start += 1
					frame_0.append(1)
				else:
					frame_0.append(0)
		return 1
	else:
		# autocorrelation numerator loop for frame
		for ii in range(len(sflag)):
			for jj in range(ii + 1, len(sflag)):
				if sflag[ii] == sflag[jj]:
					if frame_0[index] == 1:
						sigma_current += 1
				index += 1
	return sigma_current / sigma_start


def ionpair_autocorr(frame, element_0, element):
	corr = 0.0
	if len(element_0) == 0:
		return 0
	elif frame == 0:
		return 1
	else:
		for i in element:
			if i in element_0:
				corr += 1
		return corr / float(len(element_0))


def buildclusterlist(r, id2neighbors_1stshell, id2neighbors_2ndshell):  # JB
	# print "Consolidating neighbor lists into clusters..."
	autocorr = np.zeros([4, len(r)], np.float)
	avg_num = np.zeros(len(r), np.float)
	nfreeions = np.zeros(len(r), np.float)
	nSSIP = np.zeros(len(r), np.float)
	nCIP = np.zeros(len(r), np.float)
	nAGG = np.zeros(len(r), np.float)
	category_0 = [[], [], [], []]

	for t in range(len(r)):
		clusters = [[]]
		for atom_id in range(1, len(r[0])):
			# check if we've already done this atom
			already_done = False
			for cluster_id in range(1, len(clusters)):
				if atom_id in clusters[cluster_id]:
					already_done = True
					break
			if already_done:
				continue

			# if we haven't already done this atom, make it the start of a new cluster
			cluster_id = len(clusters)
			tobeadded = [atom_id]
			clusters.append(tobeadded)

			# the way this loops works is it goes through the neighbor lists adding them into the cluster, every time it finds a new atom not already in the cluster, it also checks that new atom's neighbor list for other new atoms, repeating this process until no new atoms are found
			done = False
			while not done:
				done = True
				adding = tobeadded
				tobeadded = []
				for add_id in adding:
					for neighbor_id in id2neighbors_1stshell[t][add_id]:
						if neighbor_id not in clusters[cluster_id]:
							tobeadded.append(neighbor_id)
							clusters[cluster_id].append(neighbor_id)
							done = False

		clusters.sort(key=len)

		small_clusters = [[]]
		large_clusters = [[]]
		for n in range(1, len(clusters)):
			if len(clusters[n]) <= 2:
				small_clusters.append(clusters[n])
			else:
				large_clusters.append(clusters[n])

		small_clusters.sort(key=len)
		large_clusters.sort(key=len)

		category = [[], [], [], []]
		# category: 1) free ion, 2) SSIP, 3) CIP, 4) AGG
		for atom_id in range(1, len(r[0])):
			if id2type[atom_id] == 3:  # JB
				if len(id2neighbors_2ndshell[t][atom_id]) == 0:
					nfreeions[t] += 1
					category[0].append(atom_id)
				elif len(id2neighbors_1stshell[t][atom_id]) == 0:
					nSSIP[t] += 1
					category[1].append(atom_id)
				else:
					# print t, atom_id, large_clusters
					if len(large_clusters) == 1:
						nCIP[t] += 1
						category[2].append(atom_id)
					else:
						AGG = 0
						for cluster_id in range(1, len(large_clusters)):
							if atom_id in large_clusters[cluster_id]:
								nAGG[t] += 1
								AGG = 1
								category[3].append(atom_id)
								break
						if AGG == 0:
							nCIP[t] += 1
							category[2].append(atom_id)
		# print t, category

		# reverse_clusters = [0]*len(id2type)
		# for atom_id in range(1,len(id2type)):
		#     for cluster_id in range(1,len(large_clusters)):
		#         if atom_id in large_clusters[cluster_id]:
		#             reverse_clusters[atom_id] = cluster_id
		#             break
		#     #print large_clusters

		# cluster2id.append(large_clusters)
		# id2cluster.append(reverse_clusters)

		# init sflag and num
		# sflag = []
		# for atom_id in range(1, len(r[t])):
		#     for i in range(1, len(clusters)):
		#         if atom_id in clusters[i]:
		#             sflag.append(i)    # lists each ion, characterized by a number representing the cluster they are in (in order of lowest to highest cluster number)

		# list of clusters characterized by their size (number of ions in each)
		num = np.zeros(len(clusters), np.int)
		for i in range(1, len(clusters)):
			num[i] = len(clusters[i])
		avg_num[t] = num[1:].mean()

		for count, element in enumerate(category):
			if t == 0:
				for i in element:
					category_0[count].append(i)
			autocorr[count][t] = ionpair_autocorr(t, category_0[count], element)
	return autocorr, avg_num, nfreeions, nSSIP, nCIP, nAGG


# Dump_Tools
# Modified by N. Liesen to read in unwrapped coords on 6/28/20
def read_lammpstrj_plus(fname, num_frames=float('inf'), skip_beginning=0, skip_between=0,
						flags_when_unwrap=True):  # Updated code to deal with image flags
	"""Input: fname, num_frames
    fname: filename string or 'stdin' (or a value that evaluates to false) for reading from standard in
    num_frames: optional number of frames to read before stopping, defaults to reading in all frames
    skip_beginning: skip this many frames at the beginning of the dump file
    skip_between: skip this many frames between saved frames

    Output: r, ir, timestep, box_bounds, id2type, id2mol, mol2ids
            r: num_frames by num_atoms+1 by 3 array of unscaled coordinates (indexed by frame number then atom id).
               -- Whether the coordinates are wrapped or unwrapped depends on whether wrapped or unwrapped coordinates
               are passed to the function -- NTL
            ir: num_frames by num_atoms+1 by 3 array of image flags
            timestep: num_frames length array of timesteps
            box_bounds: 3D array to store boundaries of the box, indexed by frame, x/y/z, then lower/upper
            id2type, id2mol: num_atoms+1 length arrays to map atom id to type and molecule id (if available, may be None)
            id2type[atomID], id2mol[atomID]
            mol2ids: num_mols+1 length list of atom id arrays corresponding to the molecules (if available, may be None)
                     (i.e. a list of nested 1D numpy arrays containing all atomIDs owned by the molecule) -- NTL

    NOTE: assumes that the number of atoms in the simulation is fixed
    NOTE: Only accepts wrapped and unscaled coordinates (x), wrapped and scaled coordinates (xs), or unwrapped and unscaled coordinates (xu)
          Not set up to deal with unwrapped and unscaled (xsu) coordinates or any other types. -- NTL
    NOTE: flags_when_unwrap == True returns all 0 image flags when returning unwrapped coordinates (xu) in r -- NTL
    NOTE: flags_when_unwrap == False returns None in place of ir -- NTL
    """

	def read_header(f):
		""" helper function to read in the header and return the timestep, number of atoms, and box boundaries """

		f.readline()  # ITEM: TIMESTEP
		timestep = int(float(f.readline()))

		f.readline()  # ITEM: NUMBER OF ATOMS
		num_atoms = int(float(f.readline()))

		f.readline()  # ITEM: BOX BOUNDS xx yy zz
		line = f.readline()
		line = line.split()  # xlo xhi
		xlo = float(line[0])
		xhi = float(line[1])
		line = f.readline()
		line = line.split()
		ylo = float(line[0])
		yhi = float(line[1])
		line = f.readline()
		line = line.split()
		zlo = float(line[0])
		zhi = float(line[1])

		return timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi

	# allow reading from standard input
	if not fname or fname == 'stdin':
		f = sys.stdin
	else:
		f = open(fname, 'r')

	# read in the initial header
	frame = 0
	init_timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi = read_header(f)

	# skip the beginning frames (0-->skip_beginning-1), if requested
	for skippedframe in range(skip_beginning):
		f.readline()  # ITEM: ATOMS (skip line)
		# loop over the atoms (0-->num_atoms-1) lines
		for atom in range(num_atoms):  # Skip all atom's data for selected frames
			f.readline()
		init_timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi = read_header(
			f)  # Only last read-in frame header will matter

	# preallocate arrays, if possible
	if num_frames < float('inf'):
		alloc = num_frames
		inf_frames = False
	else:
		alloc = 1
		inf_frames = True
	timestep = np.zeros(alloc, np.int)  # 1D array of timesteps
	box_bounds = np.zeros([alloc, 3, 2],
						  np.float)  # 3D array to store boundaries of the box, indexed by frame, x/y/z, then lower/upper

	timestep[frame] = init_timestep  # For frame 0, comes from header read-in -- NTL
	box_bounds[frame][0][0] = xlo
	box_bounds[frame][0][1] = xhi
	box_bounds[frame][1][0] = ylo
	box_bounds[frame][1][1] = yhi
	box_bounds[frame][2][0] = zlo
	box_bounds[frame][2][1] = zhi

	# NOTE: using num_atoms+1 here so that the arrays are indexed by their LAMMPS atom id
	r = np.zeros([alloc, num_atoms + 1, 3], np.float)  # 3D array of x, y, z coordinates, r[frame][id][coordinate]
	ir = np.zeros([alloc, num_atoms + 1, 3], np.int)  # 3D array of x, y, z image flags, r[frame][id][coordinate]

	id2mol = np.zeros(num_atoms + 1,
					  np.int)  # array to map from atom id to molecule id, builds this from the first frame, if available
	id2type = np.zeros(num_atoms + 1,
					   np.int)  # array to map from atom id to type, builds this from the first frame, if available

	# separately do the first ATOMS section so that we can initialize things, build the id2mol and id2type arrays, and so that the main loop starts with reading in the header
	line = f.readline()
	line = line.split()
	# The below determines the index for each relevant piece of info in the dump file -- NTL
	id_index = line.index(
		"id") - 2  # This whole -2 business is b/c line.split reads in ITEM: and ATOMS as the first two indices in the resulting list -- NTL
	if "mol" in line:
		mol_index = line.index("mol") - 2  # Identifies position in list containing "mol"
	else:
		mol_index = None
	if "type" in line:
		type_index = line.index("type") - 2
	else:
		type_index = None

	if "x" in line:
		scaled = False
		wrapped = True
		x_index = line.index("x") - 2
		y_index = line.index("y") - 2
		z_index = line.index("z") - 2
	elif "xs" in line:
		scaled = True
		wrapped = True
		x_index = line.index("xs") - 2
		y_index = line.index("ys") - 2
		z_index = line.index("zs") - 2
	elif "xu" in line:
		scaled = False
		wrapped = False
		x_index = line.index("xu") - 2
		y_index = line.index("yu") - 2
		z_index = line.index("zu") - 2
	else:
		print("ERROR: x coordinate not found in lammps trajectory", file=sys.stderr)
		return

	if "ix" in line:
		ix_index = line.index("ix") - 2
		iy_index = line.index("iy") - 2
		iz_index = line.index("iz") - 2
	elif "xu" in line:
		print("Coordinates are unwrapped!")
	else:
		print("ERROR: x image flag not found in lammps trajectory", file=sys.stderr)
		return

	# loop over the atoms lines for the first frame separately, the rest of the frames will be read in below
	for atom in range(num_atoms):  # num of lines in single frame output = num atoms  -- NTL
		line = f.readline()
		line = line.split()

		# get the atom id
		my_id = int(line[id_index])

		# x, y, z coordinates
		r[frame][my_id][0] = float(line[x_index])
		r[frame][my_id][1] = float(line[y_index])
		r[frame][my_id][2] = float(line[z_index])

		# unscale, if necessary
		if scaled:
			r[frame][my_id][0] = r[frame][my_id][0] * (box_bounds[frame][0][1] - box_bounds[frame][0][0]) + \
								 box_bounds[frame][0][0]
			r[frame][my_id][1] = r[frame][my_id][1] * (box_bounds[frame][1][1] - box_bounds[frame][1][0]) + \
								 box_bounds[frame][1][0]
			r[frame][my_id][2] = r[frame][my_id][2] * (box_bounds[frame][2][1] - box_bounds[frame][2][0]) + \
								 box_bounds[frame][2][0]

		# x, y, z image flags
		if wrapped:
			ir[frame][my_id][0] = int(line[ix_index])
			ir[frame][my_id][1] = int(line[iy_index])
			ir[frame][my_id][2] = int(line[iz_index])

		# if available, build the i2mol and id2type arrays
		if mol_index is not None:  # lammps may or may not print out molecule ID & type of atom ID -- NTL
			id2mol[my_id] = int(line[mol_index])
		if type_index is not None:
			id2type[my_id] = int(line[type_index])

	# build the reverse of the id2mol array
	# this is a 2D array with rows of (potentially) varying length, so nest a numpy array into a python list
	if mol_index is not None:
		num_mols = id2mol.max()  # max(molID) = total number of molecules  -- NTL
		mol2ids = [[]]  # 1-Indexed ([] pre-added) list of numpy arrays contain all atomIDs belonging to molID -- NTL
		for molid in range(1, num_mols + 1):
			# The [0] element of np.where returns an array containing the indices/atomIDs where id2mol==molid -- NTL
			mol2ids.append(np.where(id2mol == molid)[
							   0])  # Store numpy array (containing atom IDs of the selected molecule) as an entry in a list.
	else:
		num_mols = None
		mol2ids = None

	# loop over number of num_frames frames, if num_frames is infinite, will loop over all the frames in the file
	frame = 1  # this is the frame counter for frames actually read in
	frame_attempt = 0  # this is the actual frame count in the file (not counting the ones skipped in the beginning
	while frame < num_frames:

		frame_attempt += 1

		# try to read in a new header
		try:
			my_timestep, my_num_atoms, my_xlo, my_xhi, my_ylo, my_yhi, my_zlo, my_zhi = read_header(f)
		except:
			print("WARNING: hit end of file when reading in", fname, "at frame", skip_beginning + frame_attempt, file=sys.stderr)
			break

		# skip the frame if between frames to be read in and restart the loop
		if frame_attempt % (skip_between + 1) > 0:
			f.readline()  # ITEM: ATOMS
			# loop over the atoms lines
			for atom in range(num_atoms):
				f.readline()
			continue

		# if we don't know how many frames to read in, have to allocate more memory for the arrays
		if inf_frames:
			timestep = np.append(timestep, 0)

			box_bounds = np.concatenate((box_bounds, np.zeros([1, 3, 2],
															  np.float)))  # [1,3,2] --> 1 frame, 3 coords(x,y,z), 2 bounds (lower, upper) -- NTL

			r = np.concatenate((r, np.zeros([1, num_atoms + 1, 3],
											np.float)))  # [1, num_atoms+1, 3] --> 1 frame, 1-indexed atomIDs, 3 coords (x,y,z) -- NTL
			ir = np.concatenate((ir, np.zeros([1, num_atoms + 1, 3], np.float)))

		# update the timestep and box size arrays
		timestep[frame] = my_timestep
		box_bounds[frame][0][0] = my_xlo
		box_bounds[frame][0][1] = my_xhi
		box_bounds[frame][1][0] = my_ylo
		box_bounds[frame][1][1] = my_yhi
		box_bounds[frame][2][0] = my_zlo
		box_bounds[frame][2][1] = my_zhi

		f.readline()  # ITEM: ATOMS
		# loop over the atoms lines
		for atom in range(num_atoms):
			line = f.readline()
			line = line.split()

			# get the atom id
			my_id = int(line[id_index])

			# x, y, z coordinates
			r[frame][my_id][0] = float(line[x_index])
			r[frame][my_id][1] = float(line[y_index])
			r[frame][my_id][2] = float(line[z_index])

			# unscale, if necessary
			if scaled:
				r[frame][my_id][0] = r[frame][my_id][0] * (box_bounds[frame][0][1] - box_bounds[frame][0][0]) + \
									 box_bounds[frame][0][0]
				r[frame][my_id][1] = r[frame][my_id][1] * (box_bounds[frame][1][1] - box_bounds[frame][1][0]) + \
									 box_bounds[frame][1][0]
				r[frame][my_id][2] = r[frame][my_id][2] * (box_bounds[frame][2][1] - box_bounds[frame][2][0]) + \
									 box_bounds[frame][2][0]

			# x, y, z image flags
			if wrapped:
				ir[frame][my_id][0] = int(line[ix_index])
				ir[frame][my_id][1] = int(line[iy_index])
				ir[frame][my_id][2] = int(line[iz_index])

		frame += 1  # frame actually read in, increment counter -- NTL

	f.close()  # Close file

	if flags_when_unwrap and not wrapped:
		print("Flags when coordinates are unwrapped is enabled. Outputting all zero image flags.")
		return r, ir, timestep, box_bounds, id2type, id2mol, mol2ids
	elif not flags_when_unwrap and not wrapped:
		print("Flags when coordinates are unwrapped is disabled. Outputting None for image flags.")
		return r, None, timestep, box_bounds, id2type, id2mol, mol2ids
	else:
		print("Assuming coordinates are wrapped. Returning wrapped coordinates and flags.")
		return r, ir, timestep, box_bounds, id2type, id2mol, mol2ids


# End of 1st function that reads in all the data from the dump file nicely -- NTL
# r[frame, atomID, dimension]  -- contains unscaled coordiantes in units of length  (3D numpy array)
# ir[frame, atomID, dimension] -- contains periodic image flags as integers  (3D numpy array)
# box_bounds[frame, dimension, low(0)/high(1)] -- contains low/high x, y, and z bounds  (3D numpy array)
# id2type[atomID] -- corresponding integer-valued atom type for each atomID  (1D numpy array)
# id2mol[atomID] -- molecule/moleculeID for chain atom/atomID belongs to  (1D numpy array)
# timestep[frame]: -- contains timestep of frame (1D numpy array)
# mol2ids: List whose entries are numpy arrays containing all atomIDs for a selected molID. The list is indexed by molID.


# Originally from PPMD: post-pocessing polymer molecular dynamics
# Suite of tools to do post processing calculations relevant to polymers
# Originally from version for LMH's MD class fall 2016
# Original function called read_lammpstrj
# Subsuquently modified/extended by N. Liesen
# Modified N. Liesen (NTL) 2020-07-01
# Extended function to remember where it left off when you last read it. Purpose
# of this function is to be able to read a dump file in chunks for averaging coordinates
# and to be memory efficient.
def mini_read_lammpstrj(fname, num_frames=float('inf'), skip_beginning=0, skip_between=0, flags_when_unwrap=True,
						file_bookmark=None):
	def read_header(f):
		""" helper function to read in the header and return the timestep, number of atoms, and box boundaries """

		f.readline()  # ITEM: TIMESTEP
		timestep = int(float(f.readline()))

		f.readline()  # ITEM: NUMBER OF ATOMS
		num_atoms = int(float(f.readline()))

		f.readline()  # ITEM: BOX BOUNDS xx yy zz
		line = f.readline()
		line = line.split()  # xlo xhi
		xlo = float(line[0])
		xhi = float(line[1])
		line = f.readline()
		line = line.split()
		ylo = float(line[0])
		yhi = float(line[1])
		line = f.readline()
		line = line.split()
		zlo = float(line[0])
		zhi = float(line[1])

		return timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi

	# allow reading from standard input
	if not fname or fname == 'stdin':
		f = sys.stdin
	else:
		f = open(fname, 'r')

	# Resume file read from last saved point
	if not file_bookmark:
		print("No file bookmark to resume reading file from. Starting from beginning of file")
	elif isinstance(file_bookmark, int):
		print("Integer-valued file bookmark passed. Resuming file read from saved point.")
		f.seek(file_bookmark)
	else:
		sys.exit("Please pass a valid integer-valued bookmark from the tell() function.")

	# read in the initial header
	frame = 0
	init_timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi = read_header(f)

	# skip the beginning frames (0-->skip_beginning-1), if requested
	for skippedframe in range(skip_beginning):
		f.readline()  # ITEM: ATOMS (skip line)
		# loop over the atoms (0-->num_atoms-1) lines
		for atom in range(num_atoms):  # Skip all atom's data for selected frames
			f.readline()
		init_timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi = read_header(
			f)  # Only last read-in frame header will matter

	# preallocate arrays, if possible
	if num_frames < float('inf'):
		alloc = num_frames
		inf_frames = False
	else:
		alloc = 1
		inf_frames = True
	timestep = np.zeros(alloc, np.int)  # 1D array of timesteps
	box_bounds = np.zeros([alloc, 3, 2],
						  np.float)  # 3D array to store boundaries of the box, indexed by frame, x/y/z, then lower/upper

	timestep[frame] = init_timestep  # For frame 0, comes from header read-in -- NTL
	box_bounds[frame][0][0] = xlo
	box_bounds[frame][0][1] = xhi
	box_bounds[frame][1][0] = ylo
	box_bounds[frame][1][1] = yhi
	box_bounds[frame][2][0] = zlo
	box_bounds[frame][2][1] = zhi

	# NOTE: using num_atoms+1 here so that the arrays are indexed by their LAMMPS atom id
	r = np.zeros([alloc, num_atoms + 1, 3], np.float)  # 3D array of x, y, z coordinates, r[frame][id][coordinate]
	ir = np.zeros([alloc, num_atoms + 1, 3], np.int)  # 3D array of x, y, z image flags, r[frame][id][coordinate]

	# separately do the first ATOMS section so that we can initialize things, so that the main loop starts with reading in the header
	line = f.readline()
	line = line.split()
	# The below determines the index for each relevant piece of info in the dump file -- NTL
	id_index = line.index(
		"id") - 2  # This whole -2 business is b/c line.split reads in ITEM: and ATOMS as the first two indices in the resulting list -- NTL

	if "x" in line:
		scaled = False
		wrapped = True
		x_index = line.index("x") - 2
		y_index = line.index("y") - 2
		z_index = line.index("z") - 2
	elif "xs" in line:
		scaled = True
		wrapped = True
		x_index = line.index("xs") - 2
		y_index = line.index("ys") - 2
		z_index = line.index("zs") - 2
	elif "xu" in line:
		scaled = False
		wrapped = False
		x_index = line.index("xu") - 2
		y_index = line.index("yu") - 2
		z_index = line.index("zu") - 2
	else:
		print("ERROR: x coordinate not found in lammps trajectory", file=sys.stderr)
		return

	if "ix" in line:
		ix_index = line.index("ix") - 2
		iy_index = line.index("iy") - 2
		iz_index = line.index("iz") - 2
	elif "xu" in line:
		print("Coordinates are unwrapped!")
	else:
		print("ERROR: x image flag not found in lammps trajectory", file=sys.stderr)
		return

	for atom in range(num_atoms):  # num of lines in single frame output = num atoms  -- NTL
		line = f.readline()
		line = line.split()

		# get the atom id
		my_id = int(line[id_index])

		# x, y, z coordinates
		r[frame][my_id][0] = float(line[x_index])
		r[frame][my_id][1] = float(line[y_index])
		r[frame][my_id][2] = float(line[z_index])

		# unscale, if necessary
		if scaled:
			r[frame][my_id][0] = r[frame][my_id][0] * (box_bounds[frame][0][1] - box_bounds[frame][0][0]) + \
								 box_bounds[frame][0][0]
			r[frame][my_id][1] = r[frame][my_id][1] * (box_bounds[frame][1][1] - box_bounds[frame][1][0]) + \
								 box_bounds[frame][1][0]
			r[frame][my_id][2] = r[frame][my_id][2] * (box_bounds[frame][2][1] - box_bounds[frame][2][0]) + \
								 box_bounds[frame][2][0]

		# x, y, z image flags
		if wrapped:
			ir[frame][my_id][0] = int(line[ix_index])
			ir[frame][my_id][1] = int(line[iy_index])
			ir[frame][my_id][2] = int(line[iz_index])

	# loop over number of num_frames frames, if num_frames is infinite, will loop over all the frames in the file
	frame = 1  # this is the frame counter for frames actually read in
	frame_attempt = 0  # this is the actual frame count in the file (not counting the ones skipped in the beginning
	while frame < num_frames:

		frame_attempt += 1

		# try to read in a new header
		try:
			my_timestep, my_num_atoms, my_xlo, my_xhi, my_ylo, my_yhi, my_zlo, my_zhi = read_header(f)
		except:
			print("WARNING: hit end of file when reading in", fname, "at frame", skip_beginning + frame_attempt,
				  file=sys.stderr)
			break

		# skip the frame if between frames to be read in and restart the loop
		if frame_attempt % (skip_between + 1) > 0:
			f.readline()  # ITEM: ATOMS
			# loop over the atoms lines
			for atom in range(num_atoms):
				f.readline()
			continue

		# if we don't know how many frames to read in, have to allocate more memory for the arrays
		if inf_frames:
			timestep = np.append(timestep, 0)

			box_bounds = np.concatenate((box_bounds, np.zeros([1, 3, 2],
															  np.float)))  # [1,3,2] --> 1 frame, 3 coords(x,y,z), 2 bounds (lower, upper) -- NTL

			r = np.concatenate((r, np.zeros([1, num_atoms + 1, 3],
											np.float)))  # [1, num_atoms+1, 3] --> 1 frame, 1-indexed atomIDs, 3 coords (x,y,z) -- NTL
			ir = np.concatenate((ir, np.zeros([1, num_atoms + 1, 3], np.float)))

		# update the timestep and box size arrays
		timestep[frame] = my_timestep
		box_bounds[frame][0][0] = my_xlo
		box_bounds[frame][0][1] = my_xhi
		box_bounds[frame][1][0] = my_ylo
		box_bounds[frame][1][1] = my_yhi
		box_bounds[frame][2][0] = my_zlo
		box_bounds[frame][2][1] = my_zhi

		f.readline()  # ITEM: ATOMS
		# loop over the atoms lines
		for atom in range(num_atoms):
			line = f.readline()
			line = line.split()

			# get the atom id
			my_id = int(line[id_index])

			# x, y, z coordinates
			r[frame][my_id][0] = float(line[x_index])
			r[frame][my_id][1] = float(line[y_index])
			r[frame][my_id][2] = float(line[z_index])

			# unscale, if necessary
			if scaled:
				r[frame][my_id][0] = r[frame][my_id][0] * (box_bounds[frame][0][1] - box_bounds[frame][0][0]) + \
									 box_bounds[frame][0][0]
				r[frame][my_id][1] = r[frame][my_id][1] * (box_bounds[frame][1][1] - box_bounds[frame][1][0]) + \
									 box_bounds[frame][1][0]
				r[frame][my_id][2] = r[frame][my_id][2] * (box_bounds[frame][2][1] - box_bounds[frame][2][0]) + \
									 box_bounds[frame][2][0]

			# x, y, z image flags
			if wrapped:
				ir[frame][my_id][0] = int(line[ix_index])
				ir[frame][my_id][1] = int(line[iy_index])
				ir[frame][my_id][2] = int(line[iz_index])

		frame += 1  # frame actually read in, increment counter -- NTL

	file_bookmark = f.tell()  # Save file at point where we stopped reading
	f.close()  # Close file

	if flags_when_unwrap and not wrapped:
		print("Flags when coordinates are unwrapped is enabled. Outputting all zero image flags.")
		return r, ir, timestep, box_bounds, file_bookmark

	elif not flags_when_unwrap and not wrapped:
		print("Flags when coordinates are unwrapped is disabled. Outputting None for image flags.")
		return r, None, timestep, box_bounds, file_bookmark

	else:
		print("Assuming coordinates are wrapped. Returning wrapped coordinates and flags.")
		return r, ir, timestep, box_bounds, file_bookmark


# End of 1st function that reads in all the data from the dump file nicely -- NTL
# r[frame, atomID, dimension]  -- contains unscaled coordiantes in units of length  (3D numpy array)
# ir[frame, atomID, dimension] -- contains periodic image flags as integers  (3D numpy array)
# box_bounds[frame, dimension, low(0)/high(1)] -- contains low/high x, y, and z bounds  (3D numpy array)
# timestep[frame]: -- contains timestep of frame (1D numpy array)


#  Added by N. Liesen on 6/28/20 -- Adding ability to unwrap coordinates
def unwrap_coords(r, ir, box_bounds):
	"""Unwraps coordinates using wrapped coordinates and image flags. Inputs used are 3D arrays
    r[frame, atomID, axis], ir[frame, atomID, axis], and box_bounds[frame, axis, low(0)/high(1)].
    Computes xlen, ylen, and zlen (box lengths) over time using box_bounds, and unwraps coordinates.
    Function will return r_unwrap, the unwrapped coordinates as a 3D numpy array.

    In: r[frame, atomID, axis]  (3D numpy array)
        ir[frame, atomID, axis]  (3D numpy array)
        box_bounds[frame, axis, low(0)/high(1)]  (3D numpy array)
    Out: r_unwrap[frame, atomID, axis]  (3D numpy array)
    """

	box_len = box_bounds[:, :, 1] - box_bounds[:, :, 0]  # Compute box dimension time series
	xlen = box_len[:, 0]
	xlen = xlen.reshape(len(xlen), 1)  # Select x-length time series
	ylen = box_len[:, 1]
	ylen = ylen.reshape(len(xlen), 1)
	zlen = box_len[:, 2]
	zlen = zlen.reshape(len(xlen), 1)

	x_adjust = np.multiply(ir[:, :, 0], xlen)  # Broadcasting will duplicate the single column in xlen until its shape
	# matches ir[:,:,0] -- basically it will be duplicated so that there is 1 identical column per atom column in ir
	y_adjust = np.multiply(ir[:, :, 1], ylen)
	z_adjust = np.multiply(ir[:, :, 2], zlen)

	r_unwrap = np.copy(r)  # Start by coping wrapped coords, copy to avoid sharing memory
	r_unwrap[:, :, 0] = r_unwrap[:, :, 0] + x_adjust  # unwrap x-coords
	r_unwrap[:, :, 1] = r_unwrap[:, :, 1] + y_adjust
	r_unwrap[:, :, 2] = r_unwrap[:, :, 2] + z_adjust

	return r_unwrap  # Return unwrapped coordinates


# Added by N. Liesen on 6/28/20  -- required for wrap_coords
def get_box_len(box_bounds):
	""" Simple script to obtain box lengths along each axis """
	x = 0
	y = 1
	z = 2
	up = 1
	low = 0
	Lx = box_bounds[:, x, up] - box_bounds[:, x, low]  # Sigma
	Ly = box_bounds[:, y, up] - box_bounds[:, y, low]
	Lz = box_bounds[:, z, up] - box_bounds[:, z, low]
	return Lx, Ly, Lz


# Added by N. Liesen on 6/28/20  -- Adding ability to wrap coordinates
def wrap_coords(r_unwrap, bounds_box):
	"""In: r_unwrap[frame, sampleID, dimension]
           bounds_box[frame, dimension, low(0)/high(1)]
       Out: r_wrap[frame, sampleID, dimension]
            Im_flags[frame, sampleID, dimension]
    """
	x = 0
	y = 1
	z = 2
	up = 1
	low = 0

	def wrap_1D(rx, Lx, xlow):
		"""rx[sampleID],
           Lx (scalar)
           xlow (scalar)"""
		Ix = np.floor((rx - xlow) / Lx).astype(int)  # Image flag
		rx = rx - Ix * Lx  # Wrapped position
		return rx, Ix

	NUM_FRAMES = np.shape(r_unwrap)[0]
	Lx, Ly, Lz = get_box_len(bounds_box)
	xlow = bounds_box[:, x, low]
	ylow = bounds_box[:, y, low]
	zlow = bounds_box[:, z, low]

	rx = r_unwrap[:, :, x]
	ry = r_unwrap[:, :, y]
	rz = r_unwrap[:, :, z]
	r_wrap = np.zeros_like(r_unwrap, dtype=float)
	Im_flags = np.zeros_like(r_unwrap, dtype=int)

	for t in np.arange(0, NUM_FRAMES):
		r_wrap[t, :, x], Ix = wrap_1D(rx[t, :], Lx[t], xlow[t])  # Get wrapped coords & image flag
		r_wrap[t, :, y], Iy = wrap_1D(ry[t, :], Ly[t], ylow[t])
		r_wrap[t, :, z], Iz = wrap_1D(rz[t, :], Lz[t], zlow[t])

		Im_flags[t, :, x] = Ix  # For readability
		Im_flags[t, :, y] = Iy
		Im_flags[t, :, z] = Iz

	return r_wrap, Im_flags


# Added by N. Liesen on 6/28/20  -- Adding ability to scale wrapped coordinates
def scale_coords(r_wrap, bounds_box):
	"""Purpose of function is to scale wrapped coordinates.
       In: r_wrap[frame, sampleID, dimension]
           bounds_box[frame, dimension, low(0)/high(1)]
       Out: r_scale[frame, sampleID, dimension]"""
	x = 0
	y = 1
	z = 2
	up = 1
	low = 0

	def scale_1D(rx, Lx, xlow):
		"""In: rx[sampleID]
           Lx (scalar)
           xlow (scalar)"""
		rx = rx - xlow  # Shift to [0, Lx]
		return rx / Lx  # Scale to [0,1]

	NUM_FRAMES = np.shape(r_wrap)[0]
	Lx, Ly, Lz = get_box_len(bounds_box)
	xlow = bounds_box[:, x, low]
	ylow = bounds_box[:, y, low]
	zlow = bounds_box[:, z, low]

	rx = r_wrap[:, :, x]
	ry = r_wrap[:, :, y]
	rz = r_wrap[:, :, z]
	r_scale = np.zeros_like(r_wrap, dtype=float)

	for t in np.arange(0, NUM_FRAMES):
		r_scale[t, :, x] = scale_1D(rx[t, :], Lx[t], xlow[t])  # Get wrapped coords & image flag
		r_scale[t, :, y] = scale_1D(ry[t, :], Ly[t], ylow[t])
		r_scale[t, :, z] = scale_1D(rz[t, :], Lz[t], zlow[t])

	return r_scale


# Added by N. Liesen on 7/24/20  -- Adding ability to scale wrapped coordinates
# Convert xu, yu, zu --> xsu, ysu, zsu
# From LAMMPS documentation: xsu, ysu, zsu is similar to using xu, yu, zu,
# except that the unwrapped coordinates are scaled by the box size.
def scale_unwrapped_coords(r_unwrap, bounds_box, shift=True):
	""" This function scales unwrapped coordinates by the box size to
    obtain the 'xsu' style coordinates referenced in the LAMMPS dump
    documentation.
    In: r[frame, sampleID, dimension]
    bounds_box[frame, dimension, low(0)/high(1)]
    Out: r_scale[frame, sampleID, dimension]"""
	x = 0
	y = 1
	z = 2
	up = 1
	low = 0

	NUM_FRAMES = np.shape(r_unwrap)[0]
	Lx = bounds_box[:, x, up] - bounds_box[:, x, low]  # Sigma
	Ly = bounds_box[:, y, up] - bounds_box[:, y, low]
	Lz = bounds_box[:, z, up] - bounds_box[:, z, low]

	rx = r_unwrap[:, :, x]
	ry = r_unwrap[:, :, y]
	rz = r_unwrap[:, :, z]
	r_scale = np.zeros_like(r_unwrap, dtype=float)

	for t in np.arange(0, NUM_FRAMES):
		r_scale[t, :, x] = rx[t, :] / Lx[t]
		r_scale[t, :, y] = ry[t, :] / Ly[t]
		r_scale[t, :, z] = rz[t, :] / Lz[t]

	return r_scale


# Added by N. Liesen on 6/30/20  -- Req'd for correct4_center_mass
def find_beads_of_type(bead_type, id2type):
	"""Takes in the type of an atom, and the array which stores atom types (indexed by atomID)
    and generates a list of atomIDs corresponding to the selected type."""
	print("Identifying beads of type " + str(bead_type) + " and outputting relevant atomIDs.")
	atomIDs = np.where(id2type == bead_type)[0]
	return atomIDs


# Added by N. Liesen on 6/30/20  -- Req'd for correct4_center_mass
def net_mass_beads(type_bead, mass, id_to_type):
	total_mass = 0
	print("Finding beads of type " + str(type_bead))
	atomIDs = find_beads_of_type(type_bead, id_to_type)
	number_beads = len(atomIDs)
	print("number of beads of this type is " + str(number_beads))
	mass_beads = number_beads * mass
	print("net mass of these beads is " + str(mass_beads))

	return mass_beads, atomIDs


# Added by N. Liesen on 6/30/20  -- Add ability to subtract out center of mass drift from coordinates
def correct4_center_mass(r_unwrap, id_2_type, type_2_mass):
	"""In: r_unwrap[frame, atomID, dimension]"""

	num_frames = np.shape(r_unwrap)[0]

	def get_center_mass(r_t, id_to_type, type_to_mass):
		"""In: r_t[atomID, dimension]
               total_mass (scalar)
               type_to_mass (dictionary - atom types are keys)"""

		mr_t = np.zeros_like(r_t)
		total_mass = 0
		for type_bead, mass in type_to_mass.items():
			mass_beads, atomIDs_of_type = net_mass_beads(type_bead, mass, id_to_type)
			total_mass = total_mass + mass_beads
			mr_t[atomIDs_of_type, :] = r_t[atomIDs_of_type, :] * mass

		com_position_t = np.sum(mr_t, axis=0) / total_mass  # Center-of-mass position for frame t
		return com_position_t

	# Get center of mass at frame 0
	com_position_t0 = get_center_mass(r_unwrap[0, :, :], id_2_type, type_2_mass)

	# Get center of mass position for each frame and then subtract off from unwrapped coords
	r_corrected = np.zeros_like(r_unwrap)
	for t in np.arange(0, num_frames):
		com_position_t = get_center_mass(r_unwrap[t, :, :], id_2_type, type_2_mass)
		change_com = com_position_t - com_position_t0  # Find change in center-of-mass since frame 0
		# Now we want to reset the center of mass back to its frame 0 value
		r_corrected[t, :, :] = r_unwrap[t, :, :] - change_com  # Works via broadcasting - subtract change in COM

	return r_corrected  # Return unwrapped coordinates with center-of-mass reset to its frame 0 value in each frame t


# Added by N. Liesen on 7/2/20  -- Add ability to write out a lammps trajectory/dump file
def write_lammpstrj(file_name, box_bounds, timestep, id_to_mol, id_to_type, \
					r, image_flags=None, boundary_conditions=('pp', 'pp', 'pp'), \
					coordinate_type=None, write_mode='w'):
	"""Given the appropriate box boundaries, timesteps, moleculeIDs and types (for each atomID), positions,
    image flags, boundary conditions, and coordinate type/style, we can write a lammpstrajectory file, formatted
    identically to the default settings used by lammps to write trajectory files. This function doesn't return
    anything, but does write to file_name.

    In: file_name (string)
    box_bounds[frame, axis, low(0)/high(1)] (3D numpy array)
    timestep[frame] (1D numpy array)
    id_to_mol[atomID] (1D numpy array)
    id_to_type[atomID] (1D numpy array)
    r[frame, atomID, axis] (3D numpy array)
    image_flags[frame, atomID, axis] (3D numpy array)
    boundary_conditions(x, y, z) (3-tuple)
    coordinate_type (string)"""

	x = 0
	y = 1
	z = 2
	low = 0
	high = 1
	ir = image_flags

	def write_header(f, tstep, num_atoms, box_bds, coord_type):
		"""helper function to write header, timestep, number of atoms,
        and box boundaries into file"""
		xlow = box_bds[x, low]
		xhigh = box_bds[x, high]
		ylow = box_bds[y, low]
		yhigh = box_bds[y, high]
		zlow = box_bds[z, low]
		zhigh = box_bds[z, high]
		xx, yy, zz = boundary_conditions  # e.g. ('pp', 'pp', 'fs')

		f.write("ITEM: TIMESTEP\n")
		f.write('{0:0d}\n'.format(tstep))

		f.write("ITEM: NUMBER OF ATOMS\n")
		f.write('{0:0d}\n'.format(num_atoms))

		f.write('ITEM: BOX BOUNDS {} {} {}\n'.format(xx, yy, zz))
		f.write('{:.16e} {:.16e}\n'.format(xlow, xhigh))
		f.write('{:.16e} {:.16e}\n'.format(ylow, yhigh))
		f.write('{:.16e} {:.16e}\n'.format(zlow, zhigh))

		if coord_type == 'x':
			f.write('ITEM: ATOMS id mol type {} {} {} ix iy iz \n'.format('x', 'y', 'z'))
		elif coord_type == 'xs':
			f.write('ITEM: ATOMS id mol type {} {} {} ix iy iz \n'.format('xs', 'ys', 'zs'))
		elif coord_type == 'xu':
			f.write('ITEM: ATOMS id mol type {} {} {} \n'.format('xu', 'yu', 'zu'))
		elif coord_type == 'xsu':
			f.write('ITEM: ATOMS id mol type {} {} {} \n'.format('xsu', 'ysu', 'zsu'))
		else:
			sys.exit("Please input valid coordinate type: \'x\', \'xs\', \'xu\', or \'xsu\'")
		return

	def make_line(f, r_t, ir_t, atom_id, mol_id, atom_type, coord_type):
		"""Note: Adopt default lammps fomatting of %g for coordinates/floats
        and %d for integers, with all fields single space separated
        See: https://docs.python.org/2.4/lib/typesseq-strings.html
        See: https://lammps.sandia.gov/doc/dump_modify.html"""

		if coord_type == 'x' or coord_type == 'xs':
			f.write(
				"{0:0d} {1:0d} {2:0d} {3:0g} {4:0g} {5:0g} {6:0d} {7:0d} {8:0d} \n".format(atom_id, mol_id, atom_type, \
																						   r_t[x], r_t[y], r_t[z],
																						   ir_t[x], ir_t[y], ir_t[z]))
		elif coord_type == 'xu' or coord_type == 'xsu':
			f.write("{0:0d} {1:0d} {2:0d} {3:0g} {4:0g} {5:0g} \n".format(atom_id, mol_id, atom_type, \
																		  r_t[x], r_t[y], r_t[z]))
		else:
			sys.exit("Please input valid coordinate type: \'x\', \'xs\', \'xu\', or \'xsu\'")
		return

	f = file_name
	number_frames = np.shape(r)[0]
	number_atoms = np.shape(r)[1] - 1

	# allow reading from standard input
	if not file_name or file_name == 'stdin':
		f = sys.stdin
	else:
		f = open(file_name, write_mode)

	# Determine coordinate type - needs some work
	if (image_flags is None) and (coordinate_type is None):
		print("No image flags passed. Assuming unwrapped and unscaled coordinates (xu)\n")
		coordinate_type = 'xu'
	elif (not (image_flags is None)) and (coordinate_type is None):
		print("Image flags passed, but no assigned coordinate type.\n")
		print("Assuming wrapped and unscaled coordinates (x)\n")
		coordinate_type = 'x'
	elif (image_flags is None) and (coordinate_type == 'x' or coordinate_type == 'xs'):
		sys.exit("No image flags, but coordinate_type indicates wrapped coordinates.")
	elif (not (image_flags is None)) and (coordinate_type == 'xu' or coordinate_type == 'xsu'):
		print("Warning: Image flags passed, but coordinate_type indicates unwrapped coordinates.\n")
	else:  # Later logic deals with invald coordinate_type choices
		print("Coordinate type is set to {}\n".format(coordinate_type))

	print("If coordinate_type {} is incorrect, explicitly pass correct type".format(coordinate_type))
	print("(e.g. x, xs, xu, or xsu)\n")

	for t in range(0, number_frames):
		write_header(f, timestep[t], number_atoms, box_bounds[t], coordinate_type)
		for atom in range(1, number_atoms + 1):  # Write atomic coordinates
			if coordinate_type == 'x' or coordinate_type == 'xs':
				make_line(f, r[t, atom, :], ir[t, atom, :], atom, id_to_mol[atom], id_to_type[atom], coordinate_type)
			elif coordinate_type == 'xu' or coordinate_type == 'xsu':
				make_line(f, r[t, atom, :], None, atom, id_to_mol[atom], id_to_type[atom], coordinate_type)
			else:
				sys.exit("Please input valid coordinate type: \'x\', \'xs\', \'xu\', or \'xsu\'")
	f.close()
	return


# Add ability to subtract out center of mass drift from coordinates
# NOTE: By default assume x, y, and z should all be corrected
# NOTE: The below comes from the PGN bond vector field scripts
def correct4_center_mass_plus(r_unwrap, id_2_type, type_2_mass, x_is_periodic=True, y_is_periodic=True,
							  z_is_periodic=True):
	"""
    This function corrects any unwrapped coordinates which are passed for the center of mass motion which occurs over time.
    The coordinates are corrected such that the center of mass of the system matches that in the first frame in the passed
    array (r_unwrap). The function also gives the user the opportunity to decide whether the center of mass should be corrected
    along each dimenions or not using x_is_periodic, y_is_periodic, and z_is_periodic.

    Parameters
    r_unwrap[frame, atomID, dimension]: Unwrapped coordinates for all atoms in all frames; coordinates
    are indexed by the time/frame, atomID, and the x/y/z/dimension.

    id_2_type: Contains integer-valued types for each atom. The list is indexed by atomID
    (list of ints)

    type_2_mass: Contains the masses of a single bead for each atom type (dict; keys are
    integer atom types, values are float-valued masses)

    x_is_periodic, y_is_periodic, z_is_periodic: Flag indicating whether x/y/z dimension should be corrected for
    its center of mass motion. In this case whether we want to correct for the center of mass motion depends
    on the boundary conditions, but more generally these flags can be used to indicate whether the center of mass
    motion should be corrected along that dimension. (bool)

    Returns:
    r_corrected[frame, atomID, dimension]: Returns unwrapped atomic coordinates after correcting for the center of mass motion
    """

	# aliases
	x = 0
	y = 1
	z = 2

	num_frames = np.shape(r_unwrap)[0]

	def get_center_mass(r_t, id_to_type, type_to_mass):
		"""
        This function calculates the center of mass position vector for a collection of coordinates
        at a specific time/frame. This calculation requires a list of atom types for each atom, and
        a dictionary containing masses for each atom type.

        Parameters:
        r_t[atomID, dimension]: Unwrapped coordinates for all atoms at a specific time point/
        in a specific frame; indexed by frame, atomID, and the dimension/axis (3D np.array of floats)

        id_to_type: Contains integer-valued types for each atom. The list is indexed by atomID
        (list of ints)

        type_to_mass: Contains the masses of a single bead for each atom type (dict; keys are
        integer atom types, values are float-valued masses)

        Returns:
        com_position_t[dimension]: Contains the center of mass position vector for passed set
        of coordinates, r_t, after taking into account the masses of each atom type (type_to_mass)
        and the type of each atom (id_to_type) (1D np.array)
        """

		mr_t = np.zeros_like(r_t)
		total_mass = 0
		for type_bead, mass in type_to_mass.items():
			mass_beads, atomIDs_of_type = net_mass_beads(type_bead, mass, id_to_type)
			total_mass = total_mass + mass_beads
			mr_t[atomIDs_of_type, :] = r_t[atomIDs_of_type, :] * mass

		com_position_t = np.sum(mr_t, axis=0) / total_mass  # Center-of-mass position for frame t
		return com_position_t

	# Get center of mass at frame 0
	com_position_t0 = get_center_mass(r_unwrap[0, :, :], id_2_type, type_2_mass)

	# TODO: Update this function to consider whether each dimension is periodic or not
	# Get center of mass position for each frame and then subtract off from unwrapped coords
	r_corrected = np.zeros_like(r_unwrap)
	for t in np.arange(0, num_frames):
		com_position_t = get_center_mass(r_unwrap[t, :, :], id_2_type, type_2_mass)
		change_com = com_position_t - com_position_t0  # Find change in center-of-mass since frame 0
		# Now we want to reset the center of mass back to its frame 0 value
		# if x_is_periodic and y_is_periodic and z_is_periodic:
		#    r_corrected[t, :, :] = r_unwrap[t, :, :] - change_com  # works via broadcasting - subtract change in COM
		if x_is_periodic:
			r_corrected[t, :, x] = r_unwrap[t, :, x] - change_com[
				x]  # works via broadcasting - subtract change in X component of COM
		else:
			r_corrected[t, :, x] = r_unwrap[t, :, x]

		if y_is_periodic:
			r_corrected[t, :, y] = r_unwrap[t, :, y] - change_com[y]
		else:
			r_corrected[t, :, y] = r_unwrap[t, :, y]

		if z_is_periodic:
			r_corrected[t, :, z] = r_unwrap[t, :, z] - change_com[z]
		else:
			r_corrected[t, :, z] = r_unwrap[t, :, z]

	return r_corrected  # Return unwrapped coordinates with center-of-mass reset to its frame 0 value in each frame t


# Written by NLiesen during our Summer @ AFRL/RX: 6/9/21
# Provides the ability to average or smooth trajectories a la VMD
def smooth_lammpstrj(input_lammpstrj, nFrames, windowSize, windowOffset, coords_are_unwrapped, skip_frame_zero=False,
					 is_com_correction_required=False, wrap_before_writing=False, type2mass=None,
					 is_periodic=(True, True, True),
					 output_fname='smoothed.lammpstrj'):
	"""
    Determines smoothed/averaged atomic coordinates. The user is free to choose overlapping, non-overlapping,
    or moving windows for their averages through the windowSize and windowOffset settings. Smoothed atomic
    coordinates are written to a new .lammpstrj file (using our pppmd write_lammpstrj utility). Options are
    provided for the user to specify which dimensions have periodic boundaries, to choose whether to subtract
    out the center of mass motion, and to choose whether to write the averaged coordinates as wrapped or
    unwrapped. The input_lammpstrj file must contain unwrapped coordinates (TODO: Implement ability to use
    dump files with wrapped coordinates).

    NOTE: Some of the features (e.g. writing in wrapped coordinates) should be double checked before merging
    into the pppmd dump tools package.
    TODO: Add type and value checking for some of the arguments.

    e.g. Non-overlapping windows (WindowSize = 10, windowOffset = 10)
    windowID:  | --0-- | --1-- | --2-- | --3-- | --4-- | --5-- |
    frameNum: 0 1    10 11   20 21   30 31   40 41   50 51   60
                |<----->|
                windowOffset = 10
    for frame 0: skip first frame, read next 10 (frames 1-10)
    for frame N: skip first 10N + 1 frames, read next 10
    in general: skip first windowID*windowOffset + 1

    e.g. moving window average (windowSize = 10, windowOffset = 1)
    windowID: | --0-- |
             0 1      10
                | -- 1 -- |
                2         11
                  | -- 2 -- |
                  3         12
                    | -- 3 -- |
                    4         13
    for frame 0: skip first frame, read next 10 (frames 1-10)
    for frame N: skip N frames, read next 10 frames (2 - 11)
    in general: skip first windowID*windowOffset frames

    Parameters:
    input_lammptrj; .lammpstrj file containing unwrapped atomic coordinates (str)
    nFrames: Total number of frames contained in the input .lammpstrj file (TODO: make script smart enough to
    know how many frames are contained within the .lammpstjr file) (int)
    windowSize: Number of frames contained within each window/number of frames to average over for each smoothed
    frame (int)
    windowOffset: Distance between the start (or end) of window n and n-1
    coords_are_unwrapped: Does the input_lammpstrj file contain unwrapped coordinates? (bool)
    is_com_correction_required: Should the atomic coordinates be corrected for the system's COMass motion?
    skip_frame_zero: Should frame zero be omitted from the averaging process?
    wrap_before_writing: Should the outputted .lammpstrj file contain wrapped coordinates? (bool)
    type2mass: dictionary containing the masses of each atom type (dict; keys are ints, values are floats).
    Required for COM motion correction of coordinates
    is_periodic: Indicates whether the (x, y, z) coordinates are periodic or not (tuple; (bool, bool, bool))

    Returns:
    None

    However, a .lammpstrj file containing the smoothed coordinates in each non-overlapping frame is ouputted.
    """

	# handle arguments
	if type(input_lammpstrj) != str:
		raise TypeError('input lammpstrj filename must be a str')
	elif type(nFrames) != int or type(windowSize) != int or type(windowOffset) != int:
		raise TypeError('Number of frames, window size, and window offset must all be integers')

	if type(skip_frame_zero) != bool or type(is_com_correction_required) != bool \
			or type(wrap_before_writing) != bool or type(coords_are_unwrapped) != bool:
		raise TypeError('All flags must be of boolean type')

	if type(is_periodic) != tuple or len(is_periodic) != 3:
		raise TypeError('is_periodic must be passed as a tuple of length 3')

	for xyz_flag in is_periodic:
		if type(xyz_flag) != bool:
			raise TypeError('Boundary condition flags in is_periodic must be boolean')

	if is_com_correction_required and type2mass is None:
		raise Exception('Cannot correct for center of mass motion without knowing each bead type mass')

	if type2mass is not None:
		for atom_type, atom_mass in type2mass.items():
			if type(atom_type) != int:
				raise TypeError('type2mass keys must be integer-valued atomIDs')
			elif type(atom_mass) != float:
				raise TypeError('atom masses must be float-valued')

	# aliases
	timeAxis = 0
	atomAxis = 1
	dimAxis = 2
	x = 0
	y = 1
	z = 2

	# get the number of windows using the number of frames, the window size and offset
	if skip_frame_zero:
		nWindows = int(floor(((nFrames - 1.) - windowSize) / windowOffset) + 1)
	else:
		nWindows = int(floor(float(nFrames - windowSize) / windowOffset) + 1)

	print("Total number of windows: {}".format(nWindows))

	for windowID in range(nWindows):
		# read in all atomic coordinates over window of interest
		if windowID == 0:
			skip_nFrames = 0
			writeMode = 'w'
		elif windowID > 0:
			skip_nFrames = windowID * windowOffset  # skip frames to get to window of interest's start
			writeMode = 'a'  # append to previously dumped results
		if skip_frame_zero:
			skip_nFrames += 1

		print("Skipping {0} frames for window {1}".format(skip_nFrames, windowID))

		if windowID == 0:  # read_lammpstrj_plus is smart enough to detect whether coords are wrapped
			rw, _, tsteps, boxbds_window, id2type, id2mol, _ = read_lammpstrj_plus(fname=input_lammpstrj,
																				   num_frames=windowSize,
																				   skip_beginning=skip_nFrames,
																				   skip_between=0,
																				   flags_when_unwrap=False)
		elif windowID > 0:
			rw, _, tsteps, boxbds_window, _, _, _ = read_lammpstrj_plus(fname=input_lammpstrj,
																		num_frames=windowSize,
																		skip_beginning=skip_nFrames, skip_between=0,
																		flags_when_unwrap=False)

		print("Finished reading coordinates for window {}".format(windowID))

		# windowTime = (np.max(tsteps) - np.min(tsteps)) / 2.

		# average these coordinates to get the average positions of the atoms within the window
		if coords_are_unwrapped:
			if is_com_correction_required:
				# first we need to remove the center of mass motion
				rw = correct4_center_mass_plus(rw, id2type, type2mass, x_is_periodic=is_periodic[x],
											   y_is_periodic=is_periodic[y], z_is_periodic=is_periodic[z])

			rw_ave = np.mean(rw, axis=timeAxis)  # average atomic position vectors over time
			rw_ave = rw_ave.reshape((1, rw.shape[atomAxis], rw.shape[dimAxis]))
			boxbds_ave = np.mean(boxbds_window, axis=timeAxis)  # average lower and upper bounds of box over time
			boxbds_ave = boxbds_ave.reshape((1, boxbds_window.shape[1], boxbds_window.shape[2]))
			if wrap_before_writing:
				# Now we will want to re-wrap these coordinates
				rw_ave, ir_ave = wrap_coords(rw_ave, boxbds_ave)
				# Now that we have processed the required coordinates for each window we must write to a .lammpstrj
				# file.
				write_lammpstrj(output_fname, boxbds_ave,
								np.array([windowID]), id2mol, id2type, rw_ave, image_flags=ir_ave,
								coordinate_type='x', write_mode=writeMode)
			else:
				write_lammpstrj(output_fname, boxbds_ave,
								np.array([windowID]), id2mol, id2type, rw_ave, coordinate_type='xu',
								write_mode=writeMode)
		else:  # TODO: Add ability to read in wrapped coordinates using pppmd utilities
			raise Exception('Coordinates must be unwrapped prior to smoothing')

		print("Finished window number: {}".format(windowID))