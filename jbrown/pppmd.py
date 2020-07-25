#!/usr/bin/env python

# PPPMD: post-pocessing polymer molecular dynamics
# Suite of tools to do post processing calculations relevant to polymers 
# Version for LMH's MD class fall 2016
#
# Modified J. Brown 2016-10-20

import sys
import numpy as np

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
		timestep = int(f.readline())

		f.readline() # ITEM: NUMBER OF ATOMS
		num_atoms = int(f.readline())

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
		print >> sys.stderr, "ERROR: x coordinate not found in lammps trajectory"
		return

	if "ix" in line:
		ix_index = line.index("ix") - 2
		iy_index = line.index("iy") - 2
		iz_index = line.index("iz") - 2
	else:
		print >> sys.stderr, "ERROR: x image flag not found in lammps trajectory"
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
			print >> sys.stderr, "WARNING: hit end of file when reading in", fname, "at frame", skip_beginning + frame_attempt
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
	for t in xrange(frames):
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
	for t in xrange(frames):
		# loop over all atoms
		for atom_id in xrange(1,num_atoms):
			# now find the distance from this atom to all the following atoms 
			# using numpy to vectorize the calculation for speed
			delta_r = np.abs(r[t][atom_id+1:] - r[t][atom_id]) # creates an array of vectors containing the absolute differences in the x, y, and z directions
			delta_r = np.where(delta_r > 0.5 * box_size, delta_r - box_size, delta_r) # this does an element by element comparison, wherever the absolute difference is greater than half the box size, we should have used a periodic image of one of the two atoms, so to correct this, shift that value by the appropriate box dimension. Note that since this happens on an element by element basis, it may shift the difference in the x-direction by the size of the box in the x-direction, but leave the rest of the vector alone.
			dist = np.sqrt((delta_r ** 2).sum(axis=1)) # now calculate the distances, the axis=1 option makes the sums happen over the vectors not the entire dataset, so the resulting data is an array containg the distances from atom_id
			
			# bin the distances into the gofr array
			dist = dist/bin_size
			for n in xrange(len(dist)):
				bin_num = int(np.floor(dist[n]))
				if bin_num < num_bins:
					gofr[bin_num] += 1.0

	# rescale based on number of atoms, bin volume, density, and number of frames
	for bin_num in xrange(num_bins):
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
	for t in xrange(frames):
		# loop over a grid of k-vectors inside the sphere of all indices that makes all k-vectors within the maximun range
		for ix in xrange(num_bins + 1):
			k_x = ix*dk[0]
			for iy in xrange( int(np.floor(np.sqrt(num_bins**2 - ix**2))) + 1):
				k_y = iy*dk[1]
				for iz in xrange( int(np.floor(np.sqrt(num_bins**2 - ix**2 - iy**2))) + 1):
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
