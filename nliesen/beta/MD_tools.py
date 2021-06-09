#!/usr/bin/env python2
# coding: utf-8

# Purpose: Original focus on getting neighbors within some spherical cutoff of each atom
# Author: N. Liesen(NTL)
# Resources: Most of the code is stolen from cluster_v3.1.py, which was a collaborative
# script from the Hall group (Jon/Anne/Kyaw/Janani/Conner and others?) and included the
# main code required for neighbor listing
# Date Created: 6/27/20 -- NTL
# Last Modified: 8/20/20 -- NTL

import numpy as np
import random
#import os
#import sys
import matplotlib.style
import matplotlib.font_manager
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
plt.interactive(True)
from mpl_toolkits.mplot3d import Axes3D
from math import pi, sqrt, cos, sin, floor

# aliases
x = 0
y = 1
z = 2
low = 0
high = 1

#%%

# Functions for generating test data

# required for wrap_coords
# from pppmd nliesen subpackage
def get_box_len(box_bounds):
    """ Simple script to obtain box lengths along each axis """
    x = 0
    y = 1
    z = 2
    up = 1
    low = 0
    Lx = box_bounds[:, x, up]-box_bounds[:, x, low]  # Sigma
    Ly = box_bounds[:, y, up]-box_bounds[:, y, low]
    Lz = box_bounds[:, z, up]-box_bounds[:, z, low]
    return Lx, Ly, Lz


#from pppmd nliesen subpackage
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
        Ix = np.floor((rx-xlow)/Lx).astype(int)  # Image flag
        rx = rx - Ix*Lx  # Wrapped position
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


def make_random_vector(desired_radius = 1.0):
    """This script will generate a vector (vx, vy, vz) which points in a
    random direction on the surface of a sphere of the specified radius. 
    By default this will generate a vector on the surface of a unit sphere.
    refs: https://mathworld.wolfram.com/SpherePointPicking.html
    https://math.stackexchange.com/questions/44689/how-to-find-a-random-axis-or-unit-vector-in-3d
    """
    theta = random.random()*2*pi  # theta dist uniform over [0, 2pi]
    vz = random.random()*2 - 1  # uniform dist over cos(phi) from [-1,1]
    vy = sqrt(1-vz**2)*sin(theta)
    r = sqrt(vx**2+vy**2+vz**2)
    scale = desired_radius/r  # Normalize vector
    return scale*np.array([vx, vy, vz])  # return unit vector

#%%

# Neighbor listing functions

def get_box_length(box_bounds, frame):
    x = 0
    y = 1
    z = 2
    low = 0
    high = 1
    # break box into bins
    Lx = box_bounds[frame][x][high] - box_bounds[frame][x][low]
    Ly = box_bounds[frame][y][high] - box_bounds[frame][y][low]
    Lz = box_bounds[frame][z][high] - box_bounds[frame][z][low]
    return Lx, Ly, Lz


def get_number_bins(L, RCUT):
    """Determines smallest use-able dimensions for cells in neighbor listing process, using
     box length, L, and neighbor list cutoff, RCUT."""
    # Find number of bins
    number_bins = int(floor(L / RCUT))  # floor b/c L > RCUT, will only check adjacent bins
    size_bin = L / float(number_bins)
    return number_bins, size_bin


def print_bin_info(num_bins, size_bin):
    print "{0:0d} bins along the x-axis of size {1:0g}".format(num_bins, size_bin)
    return


def is_on_same_molecule(atomID_a, atomID_b, id_to_mol):
    moleculeID_a = id_to_mol[atomID_a]
    moleculeID_b = id_to_mol[atomID_b]
    return moleculeID_a == moleculeID_b

def atomIDs_too_close(atomID_a, atomID_b, threshold):
    # For use if atoms are on the same chain, and IDs are sequential
    beads_between = np.abs(atomID_a - atomID_b) - 1
    return beads_between < threshold


def get_molecule_lens(molecule_to_ids, more_beads_then=2):
    """Find length of grafted chains using a list containing a list of atomIDs belonging
    to each molID, molecule_to_ids. This is done by finding the atom-based lengths associated
    with each moleculeID and then excluding molecules which are too short (more_beads_than)
    and then checking if all entries are of equal length. In the monodisperse case we get
    the uniform lengths of the grafted chains.
    Input: molecule_to_ids: List of 
    """
    if len(molecule_to_ids) < 1:  # Make sure list passed isn't empty
        raise Exception('molecule_to_ids list is empty')

    all_molecule_lens = [len(molecule) for molecule in molecule_to_ids]
    molecule_lens = [len(molecule) for mol in all_molecule_lens if mol > more_beads_then]
    if len(set(molecule_lens)) == 1:  # Check if all entries are the same length
        molecule_length = molecule_lens[0]  # Pull out chain lengths
        return molecule_length
    else:
        raise Exception('Are your grafted chains polydisperse/different lengths?')


def get_bin_dimensions(Lx, Ly, Lz, cutoff, sel_frame):
    number_xbins, size_xbin = get_number_bins(Lx, cutoff)
    number_ybins, size_ybin = get_number_bins(Ly, cutoff)
    number_zbins, size_zbin = get_number_bins(Lz, cutoff)

    print_bin_info(number_xbins, size_xbin)
    print_bin_info(number_ybins, size_ybin)
    print_bin_info(number_zbins, size_zbin)

    return number_xbins, number_ybins, number_zbins, size_xbin, size_ybin, size_zbin

# Requires: get_box_length
def bin_atoms_in_3D(number_bins, box_bound, positions, chosen_frame):
    # FIXME: When getting the number of atoms we are accidently using the global
    # r_wrap
    # FIXME: Accidently using global box_bounds instead of box_bound when
    # getting box lengths
    """This function will break the simulation domain into 3D cells ("voxels") using the
    specified number_bins = (num_xbins, num_ybins, num_zbins) list and then
    assign each atom to a voxel depending on the position coords, positions of rr, and
    the box dimensions, box_bound, passed to the function. We also need to pass a frame
    index since box_bound = box_bound[frame, axis, low(0)/high(1)] and
    rr = rr[frame, atomID, axis]. Once each atom is binned for the chosen frame
    the function will return a 3D list bin_to_id indexed by the x, y, and z bin indices,
    and a 1D list of (nested) lists, id_to_bin, containing the x/y/z or i/j/k indices of
    the bin to which the selected atomID was assigned.
    Input: positions[frame, atomID, axis]: contains all atom position vectors over time (3D np array)
           box_bound[frame, axis, low(0)/high(1)]: contains max/min box bounds for x/y/z axis over time (3D np array)
           number_bins[axis]: contains number of bins along each axis (1D list)
           chosen_frame: int-valued frame at which to do analysis (int)
    Output: bin_to_id[xbin, ybin, zbin]: contains atom IDs assigned to each x, y, and z-index specified bin (3D list)
            id_to_bin[atomID]: contains a nested list w/ x, y, and z-index for cell it is assigned to (1D list) 
    """
    #aliases
    x = 0
    y = 1
    z = 2
    low = 0
    high = 1
    rr = positions
    tt = chosen_frame

    #number_atoms = np.shape(r_wrap)[1] - 1
    number_atoms = np.shape(rr)[1] - 1
    num_xbins = number_bins[x]
    num_ybins = number_bins[y]
    num_zbins = number_bins[z]
    bin_to_id = [[[ [] for k in xrange(num_zbins)] for j in xrange(num_ybins)] for i in xrange(num_xbins)]  # for atomIDs owned by voxel
    id_to_bin = [ [] for n in xrange(number_atoms + 1) ]  # for each atomID - tells us which voxel it's in

    # needed to shift xlow --> xhigh (where xhigh = xlow+Lx) to 0 --> Lx
    xlow = box_bound[tt][x][low]
    ylow = box_bound[tt][y][low]
    zlow = box_bound[tt][z][low]

    # map the atom locations into those 3D bins
    #Lx, Ly, Lz = get_box_length(box_bounds, tt)
    Lx, Ly, Lz = get_box_length(box_bound, tt)
    for atom_id in xrange(1, number_atoms + 1):
        #TODO: should I keep this 1e-10? - was in original clustering/nlist code from Dr. Hall.
        #NOTE: "It's easy to keep it and make it smaller, and see if anyhing changes"
        i = int(floor( num_xbins*( (rr[tt][atom_id][x] - xlow) / (Lx+1e-10) ) ))
        j = int(floor( num_ybins*( (rr[tt][atom_id][y] - ylow) / (Ly+1e-10) ) ))
        k = int(floor( num_zbins*( (rr[tt][atom_id][z] - zlow) / (Lz+1e-10) ) ))
        #i = int(floor( num_xbins*( (rr[tt][atom_id][x] - xlow) / Lx ) ))
        #j = int(floor( num_ybins*( (rr[tt][atom_id][y] - ylow) / Ly ) ))
        #k = int(floor( num_zbins*( (rr[tt][atom_id][z] - zlow) / Lz ) ))
        bin_to_id[i][j][k].append(atom_id)
        id_to_bin[atom_id] = [i, j, k]

    return bin_to_id, id_to_bin


class Neighbors:
    """Each neighbor will have it's atomID stored in self.atom_ids. The atomID's corresponding molID will be stored in
    self.mol_ids. The distance used to determine that the atom belongs in this neighbor list is stored in self.distances.
    The add_neighbor method will be called to add the atomID, molID, and distance into the proper list attribute."""
    def __init__(self,  central_atomID):
        self.center_atomID = central_atomID  # To verify the neighbor list is 4 correct atom
        self.atom_ids = []
        self.mol_ids = []
        self.distances = []

    def add_neighbor(self, new_atomID, new_molID, dist):
        self.atom_ids.append(new_atomID)
        self.mol_ids.append(new_molID)
        self.distances.append(dist)
        return

    def __str__(self):  # modify print behavior of class instances
        number_neighbors = len(self.atom_ids)
        nneigh_dist = np.min(np.array(self.distances))
        farthest_dist = np.max(np.array(self.distances))
        return "atomID {0} has {1} neighbors. ".format(self.center_atomID, number_neighbors) + \
            "The nearest and farthest neighbors are at a distance of {0} and {1} respectively.".format(nneigh_dist, farthest_dist)


#%%

# Req'd functions: get_box_length, get_bin_dimensions, bin_atoms_in_3D
# Req'd classes: Neighbors
# Note: The special_neighs feature relies on atomIDs being sequential in the traj file!! 
def make_neighbor_list(r_wrap, box_bounds, RCUT, frame, atom2molid, special_neighs=3):
    # FIXME: id2mol should be an argument, right now it is grabbing a global variable
    """This function will break the simulation domain into 3D cells ("voxels") and then
    assign each atom to a voxel depending on the wrapped position coordinates, r_wrap,
    passed to the function. The size of the side length of cells is decided based on
    the radial cutoff distance, RCUT, used when assigning atoms to neighbor lists. The
    idea is to set the cell size such that when looking for all neighbors within RCUT
    we need only look into the 26 voxels surrounding the central voxel (in which the
    central atom is contained). This is a cheap way to minimize the effort required to
    build the neighbor list. The minimum possible distance is to set the side length
    of the cells to be greater than 1/2(RCUT). The last step of the algorithm is to
    check distances between pairs of atoms, in order to decide whether they are neighbors.
    
    Input: r_wrap[frame, atomID, axis]: atomic position vectors over time for all atoms (3D np array)
    box_bounds[frame, axis, low(0)/high(1)]: time series containing min/max box bounds for each axis (3D np array) 
    RCUT: radial distance below which atoms are included in the central atoms neighbor list (float)
    frame: Specifies which snapshot to analyze (int)
    special_neighs: Number of beads between 2 beads on the same chain before we count the interactions
                    For instance, 1-2, 1-3, and 1-4 bead interactions may not be included in the neighbor list
                    You can also use special_neighs=chain_length-1 if you want to completely exclude any beads
                    on the same chain from the list (int: default 3)
    
    Output: beads: List whose entries are neighbor objects for each atom in the system. The neighbor object
                   contains information about the atom being considered, all its neighbors, the molecule IDs of
                   each of its neighbors, and the distance between the bead of interest and each neighbor.
                   id2neighbors[some_atomID]: 1D list of Neighbor objects containing all information about each atom's neigh list 
                                    id2neighbors[some_atomID].atom_ids: atom_ids in some_atomID's neigh list (1D list)
                                    id2neighbors[some_atomID].distances: distances for each atom ID in neigh list (1D list)
                                    id2neighbors[some_atomID].mol_ids: molecule ids for each atom ID in neigh list (1D list)
    """

    number_atoms = np.shape(r_wrap)[1] - 1
    Lx, Ly, Lz = get_box_length(box_bounds, frame)
    num_xbins, num_ybins, num_zbins, _, _, _ = get_bin_dimensions(Lx, Ly, Lz, RCUT, frame)
    print "Started 3D binning procedure..."
    num_bins = [num_xbins, num_ybins, num_zbins]
    bin2id, id2bin = bin_atoms_in_3D(num_bins, box_bounds, r_wrap, frame)
    print "Done binning!"

    print "Building neighbor list..."

    id2neighbors = [ Neighbors(atom_i) for atom_i in range(0, number_atoms + 1)]  # List of neighbor objects
    id2neighbors[0] = None  # atomIDs are 1-indexed

    # loop over beads and their local area to build a neighbor list
    for atom_i in xrange(1, number_atoms + 1):

        atoms_xbin = id2bin[atom_i][x]
        atoms_ybin = id2bin[atom_i][y]
        atoms_zbin = id2bin[atom_i][z]
        bins_indices = [str(atoms_xbin), str(atoms_ybin), str(atoms_zbin)]

        if atom_i%10000 == 0:
            print "Getting neighbors for atom with id " + str(atom_i)
            print "Looping over the 26 bins surrounding bin of index " + ', '.join(bins_indices) + " and the bin itself."

        # loop over the local area of each bead (all bins touching the atom's bin), and shift at the periodic BCs
        first_bin = 0
        last_xbin = num_xbins - 1
        last_ybin = num_ybins - 1
        last_zbin = num_zbins - 1

        # Loops over 27 cells/voxels - central(atomID belongs to this) + touching cells
        for i in xrange(atoms_xbin-1, atoms_xbin+2):
            if i == -1:
                i = last_xbin  # Periodic
                xshift = -Lx
            elif i == last_xbin + 1:
                i = first_bin  # Periodic
                xshift = Lx
            else:
                xshift = 0

            for j in xrange(atoms_ybin-1, atoms_ybin+2):
                if j == -1:
                    j = last_ybin
                    yshift = -Ly
                elif j == last_ybin + 1:
                    j = first_bin
                    yshift = Ly
                else:
                    yshift = 0

                for k in xrange(atoms_zbin-1, atoms_zbin+2):
                    if k == -1:
                        k = last_zbin
                        zshift = -Lz
                    elif k == last_zbin + 1:
                        k = first_bin
                        zshift = Lz
                    else:
                        zshift = 0

                    # loop over the beads in this box and calculate the distance
                    for test_id_j in bin2id[i][j][k]:
                        # if it's already in the nhood, don't bother
                        if test_id_j in id2neighbors[atom_i].atom_ids:  # check if test atom is already in atom i's nlist
                            continue
                        if test_id_j == atom_i:  # Check that the two atoms aren't the same
                            continue
                        if is_on_same_molecule(atom_i, test_id_j, atom2molid):  # Are atoms on the same molecule?
                            if atomIDs_too_close(atom_i, test_id_j, special_neighs):  # Are there at least special_neighs beads between the 2 atoms
                                continue

                        dx = r_wrap[frame][test_id_j][x] + xshift - r_wrap[frame][atom_i][x]
                        dy = r_wrap[frame][test_id_j][y] + yshift - r_wrap[frame][atom_i][y]
                        dz = r_wrap[frame][test_id_j][z] + zshift - r_wrap[frame][atom_i][z]
                        dr2 = dx*dx + dy*dy + dz*dz

                        # if close enough add to the neighbor list
                        if dr2 < RCUT**2:
                            id2neighbors[atom_i].add_neighbor(test_id_j, atom2molid[test_id_j], np.sqrt(dr2))
                            id2neighbors[test_id_j].add_neighbor(atom_i, atom2molid[atom_i], np.sqrt(dr2))
    return id2neighbors  #return instance of Atoms class containing list of Neighbor objects w/ atomIDs, molIDs, and distances for neighbors


# %%

if __name__ == '__main__':  # set up to run as jupyter notebook cells
    # This test code will generate X atoms on the surface of spheres of several different radii around a
    # central atom (my_center), will then generate a neighbor list and plot all the atoms, drawing red
    # Xs over atoms which are contained within the neighbor list of the central atom. The resulting plot
    # is 3D and interactive provided appropriate packages are installed for widget mode in jupyter.
    % matplotlib widget
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    #for r_atom in r_wrap[0]:
    #    ax.scatter3D(r_atom[x], r_atom[y], r_atom[z])

    # Pick an atom and generate position vectors at a specific radius
    # away from the sphere
    my_center = np.array([5, 5, 5])
    ax.scatter3D(my_center[x], my_center[y], my_center[z])

    nFrames = 1
    my_radii = [2.0, 3.99, 10]
    nRadii = len(my_radii)
    points_per_radius = 100
    nAtoms = nRadii * points_per_radius + 1  # + 1 for the center atom
    nDims = 3
    r = np.zeros((nFrames, nAtoms + 1, nDims))
    id2mol = [0]

    # generate atomic positions such that the atoms fall onto the surface
    # of spheres of different radii
    frame = 0
    atomid = 1 
    r[frame, atomid] = my_center
    id2mol.append(atomid)
    atomid += 1
    for radius in my_radii:
        for ii in range(points_per_radius):
            my_vector = make_random_vector(desired_radius = radius)
            r[frame, atomid] = my_center + my_vector 
            id2mol.append(atomid)  # new molecule for every placed bead
            atomid += 1

    # make box bounds for data

    # test w/out worrying about PBCs
    xmin = my_center[x] - np.max(np.array(my_radii))
    xmax = my_center[x] + np.max(np.array(my_radii))

    # test when wrapped through box edge
    #xmin = my_center[x]
    #xmax = my_center[x] + 2*np.max(np.array(my_radii))

    # test w/out worrying about PBCs
    ymin = my_center[y] - np.max(np.array(my_radii))
    ymax = my_center[y] + np.max(np.array(my_radii))

    # test when wrapped through box edge
    #ymin = my_center[y]
    #ymax = my_center[y] + 2*np.max(np.array(my_radii))

    # test w/out worrying about PBCs
    zmin = my_center[z] - np.max(np.array(my_radii))
    zmax = my_center[z] + np.max(np.array(my_radii))

    # test when wrapped through box edge
    #zmin = my_center[z]
    #zmax = my_center[z] + 2*np.max(np.array(my_radii))

    boxbds = np.zeros((nFrames, nDims, 2))
    boxbds[frame, x, low] = xmin
    boxbds[frame, x, high] = xmax
    boxbds[frame, y, low] = ymin
    boxbds[frame, y, high] = ymax
    boxbds[frame, z, low] = zmin
    boxbds[frame, z, high] = zmax

    # wrap the position vectors
    r, ir = wrap_coords(r, boxbds)

    # Plot the generated atom positions
    for atom_positions in r[frame, :, :]:
        ax.scatter3D(atom_positions[x], atom_positions[y], atom_positions[z])

    # Now that we have atomic positions we can test the neighbor list
    id2neighlist = make_neighbor_list(r, boxbds, 4, frame, id2mol, special_neighs=0)

    neighbors = id2neighlist[1].atom_ids
    neighbors.sort() 
    print neighbors

    for beadID in neighbors:
        ax.scatter3D(r[frame, beadID, x], r[frame, beadID, y], r[frame, beadID, z], marker = 'x', color='red')

    print(id2neighlist[1])

    # Check that if atom i's nlist contains j, that j's nlist contains i
    for atom_i in range(1, nAtoms + 1):
        for atom_j in id2neighlist[atom_i].atom_ids:
            if atom_i not in id2neighlist[atom_j].atom_ids:
                print ("Oh no!")

    print "Done"
