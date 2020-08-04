#!/usr/bin/python
 
# Script:  residencetime.py
# Purpose: calculate ion pair correlation function
# Syntax:  python residencetime.py 
# Author:  Kevin Shen, May 2019

import sys
import numpy as np 
import pandas as pd

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
    print "Reading configurations..."
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
        print "Skipping " + str(skippedframe) 
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

    # NOTE: using num_atoms+1 here so that the arrays are indexed by their LAMMPS atom id
    r_temp = np.zeros([num_atoms+1, 3], np.float) # 3D array of x, y, z coordinates, r[frame][id][coordinate]
    ir_temp = np.zeros([num_atoms+1, 3], np.int) # 3D array of x, y, z image flags, r[frame][id][coordinate]

    # id2mol = np.zeros(num_atoms+1, np.int) # array to map from atom id to molecule id, builds this from the first frame, if available
    index2type = np.zeros(num_atoms+1, np.int) # array to map from atom id to type, builds this from the first frame, if available
    index2id = np.zeros(num_atoms+1, np.int) # array to map from atom id to type, builds this from the first frame, if available
    
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
            r_temp[my_id][0] = r_temp[my_id][0]*(box_bounds[frame][0][1]-box_bounds[frame][0][0]) + box_bounds[frame][0][0]
            r_temp[my_id][1] = r_temp[my_id][1]*(box_bounds[frame][1][1]-box_bounds[frame][1][0]) + box_bounds[frame][1][0]
            r_temp[my_id][2] = r_temp[my_id][2]*(box_bounds[frame][2][1]-box_bounds[frame][2][0]) + box_bounds[frame][2][0]

        # x, y, z image flags
        ir_temp[my_id][0] = int(line[ix_index])
        ir_temp[my_id][1] = int(line[iy_index])
        ir_temp[my_id][2] = int(line[iz_index])

    # NOTE: using num_types+1 here so that the arrays are indexed by their LAMMPS atom id
    num_not_types = num_atoms - num_types
    r = np.zeros([alloc, num_types+1, 3], np.float) # 3D array of x, y, z coordinates, r[frame][id][coordinate]
    ir = np.zeros([alloc, num_types+1, 3], np.int) # 3D array of x, y, z image flags, r[frame][id][coordinate]

    id2type = np.zeros(num_types+1, np.int) # array to map from atom id to type, builds this from the first frame, if available
    id2index =  np.zeros(num_types+1, np.int) # array to map from atom id to index, builds this from the first frame, if available

    # store the temporary data into real arrays
    my_id = 0
    for atom in range(num_atoms):
        index = atom+1
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

            r = np.concatenate( ( r, np.zeros([1, num_types+1, 3], np.float) ) )
            ir = np.concatenate( ( ir, np.zeros([1, num_types+1, 3], np.float) ) )
        
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
            index = int(line[id_index])
            if index2type[index] in types:
                my_id = index2id[index]

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
    
        print "Reading frame {}".format(frame)
        frame += 1
    
    print '===============Summary==============='
    print 'Total number of atoms =', num_atoms
    print 'Total number of ions =', num_types

    return r, ir, timestep, box_bounds, id2type, id2index


def buildnlist(r, bin2id, id2bin, bins, boxsize, id2type, dist_range, nearest):
    print 'Building neighbor lists...'
    #loop over the local area of each bead, and shift at the periodic BCs
    id2neighbors = [[[] for n in xrange(len(r[0]))]for n in xrange(len(r))]
    for t in range(len(r)):
        # loop over atoms
        for atom in range(1, len(r[0])):
            shift = np.zeros(3, np.float)
            binloc = id2bin[t][atom]
            if nearest:
                dmin = dist_range**2
                dmin_id = 0
            for i in range(binloc[0]-1,binloc[0]+2):
                if i == -1:
                    i = bins[0]-1
                    shift[0] = -boxsize[0]
                elif i == bins[0]:
                    i = 0
                    shift[0] = boxsize[0]
                else:
                	shift[0] = 0
            
                for j in range(binloc[1]-1,binloc[1]+2):
                    if j == -1:
                        j = bins[1] - 1
                        shift[1] = -boxsize[1]
                    elif j == bins[1]:
                        j = 0
                        shift[1] = boxsize[1]
                    else:
                		shift[1] = 0

                    for k in range(binloc[2]-1,binloc[2]+2):
                        if k == -1:
                            k = bins[2] - 1
                            shift[2] = -boxsize[2]
                        elif k == bins[2]:
                            k = 0
                            shift[2] = boxsize[2]
                        else:
                			shift[2] = 0
                            
                        #loop over the beads in this box and calculate distance
                        for test_id in bin2id[t][i][j][k]:
                            if (test_id in id2neighbors[t][atom]) or (test_id == atom) or (id2type[test_id] == id2type[atom]):
                                continue         
                            dr = r[t][test_id] + shift - r[t][atom]
                            dr2 = dr.dot(dr)
                            if dr2 < dist_range**2:
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
def binbox(r, boxsize, bound_lo, dist_range):
    print 'Binning box...'
    bins = np.floor(boxsize/dist_range).astype(int)
    id2bin = np.zeros([len(r), len(r[0]), 3], np.int)
    bin2id = [[[[[] for k in range(bins[2])] for j in range(bins[1])] for i in range(bins[0])] for t in range(len(r))]

    # loop over frames
    for t in range(len(r)):
        # loop over atoms
        for atom in range(1, len(r[0])):
            xbin, ybin, zbin = np.floor((r[t][atom] - bound_lo)/boxsize*bins).astype(int)
            try:
                id2bin[t][atom] = (xbin, ybin, zbin)
                bin2id[t][xbin][ybin][zbin].append(atom)
            except:
                print >> sys.stderr, 'Warning: atom {} at timestep {} slightly outside of box: '.format(atom, t)
                print >> sys.stderr, '{} > {}'.format(r[t][atom], bound_hi)

    return bin2id, id2bin, bins

# NOTE: Based on https://lammps.sandia.gov/threads/msg32219.html, atoms can be slightly outside the periodic boundary
# if such thing happens, we wrap it back
def wrap(r, box_bounds):
    print 'Wrapping...'
    bound_lo = np.array([ box_bounds[0][0][0], box_bounds[0][1][0], box_bounds[0][2][0] ])
    bound_hi = np.array([ box_bounds[0][0][1], box_bounds[0][1][1], box_bounds[0][2][1] ])
    boxsize = bound_hi-bound_lo
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
                elif id2neighbors[t][atom] == id2neighbors[0][atom] and corr[t-1][atom] != 0:
                    corr[t][atom] = 1
    return corr.mean(axis=1)/h.mean()


def cluster_autocorr(frame, sflag):
    global frame_0, sigma_start
    sigma_current = 0.0  #numerator of ACF
    index = 0

    # establish where ions are at the beginning
    if frame == 0:
        frame_0 = []
        sigma_start = 0.0    #denominator of ACF
        for ii in range(len(sflag)):
            for jj in range(ii+1, len(sflag)):
                if sflag[ii] == sflag[jj]:
                    sigma_start += 1
                    frame_0.append(1)
                else:
                    frame_0.append(0)
        return 1
    else:
        # autocorrelation numerator loop for frame
        for ii in range(len(sflag)):
            for jj in range(ii+1, len(sflag)):
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
        return corr/float(len(element_0))

def buildclusterlist(r, id2neighbors_1stshell, id2neighbors_2ndshell):
    # print "Consolidating neighbor lists into clusters..."
    autocorr = np.zeros([4, len(r)], np.float)
    avg_num = np.zeros(len(r), np.float)
    nfreeions = np.zeros(len(r), np.float)
    nSSIP = np.zeros(len(r), np.float)
    nCIP = np.zeros(len(r), np.float)
    nAGG = np.zeros(len(r), np.float)
    category_0 = [[],[],[],[]]

    for t in range(len(r)):
        clusters  = [ [] ]
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
                done  = True
                adding = tobeadded
                tobeadded = []
                for add_id in adding:
                    for neighbor_id in id2neighbors_1stshell[t][add_id]:
                        if neighbor_id not in clusters[cluster_id]:
                            tobeadded.append(neighbor_id)
                            clusters[cluster_id].append(neighbor_id)
                            done = False
            
        clusters.sort(key=len)
        
        small_clusters = [ [] ]
        large_clusters = [ [] ]
        for n in range(1, len(clusters)):
            if len(clusters[n]) <= 2:
                small_clusters.append(clusters[n])
            else:
                large_clusters.append(clusters[n])

        small_clusters.sort(key=len)
        large_clusters.sort(key=len)

        category = [[],[],[],[]]
        # category: 1) free ion, 2) SSIP, 3) CIP, 4) AGG
        for atom_id in range(1, len(r[0])):
            if id2type[atom_id] == 3:
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
                        for cluster_id in range(1,len(large_clusters)):
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
        #     for cluster_id in xrange(1,len(large_clusters)):
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

