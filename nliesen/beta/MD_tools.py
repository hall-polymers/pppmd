#!/usr/bin/env python
# coding: utf-8

# Purpose: Original focus on getting neighbors within some spherical cutoff of each atom
# Author: N. Liesen(NTL)
# Date Created: 6/27/20 -- NTL
# Last Modified: 8/20/20 -- NTL

import numpy as np
from math import floor, sqrt
import sys 
import os


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
        sys.exit('molecule_to_ids list is empty')
    
    all_molecule_lens = [len(molecule) for molecule in molecule_to_ids]
    molecule_lens = [len(molecule) for mol in all_molecule_lens if mol > more_beads_then]
    if len(set(molecule_lens)) == 1:  # Check if all entries are the same length
        molecule_length = molecule_lens[0]  # Pull out chain lengths
        return molecule_length
    else:
        sys.exit('Are your grafted chains polydisperse/different lengths?')


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
    
    number_atoms = np.shape(r_wrap)[1] - 1
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
    Lx, Ly, Lz = get_box_length(box_bounds, tt)
    for atom_id in xrange(1, number_atoms + 1):  #TODO: should I keep this 1e-10? - was in original clustering/nlist code from Dr. Hall.    
        i = int(floor( num_xbins*( (rr[tt][atom_id][x] - xlow) / (Lx+1e-10) ) ))
        j = int(floor( num_ybins*( (rr[tt][atom_id][y] - ylow) / (Ly+1e-10) ) ))
        k = int(floor( num_zbins*( (rr[tt][atom_id][z] - zlow) / (Lz+1e-10) ) ))
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


class Atoms:  # Store neighbor objects with atom_IDs, mol_IDs, and distances
    """This object will contain a list (indexed by atomID) containing a Neighbor type object. These stored
    objects will contain all relevant information about the chosen neighbor."""
    def __init__(self, num_atoms):
        self.id2neighs = [ None for i in range(0, num_atoms + 1)]  # Neighbors objects indexed by atomID
        
    def store_neighs_of(self, central_atomID, neighbor_list):
        self.id2neighs[central_atomID] = neighbor_list  # Store nlist obj in atom2neighs
        return
                            
    def get_neighs_of(self, central_atomID):
        return self.id2neighs[chosen_atomsID]


# Req'd functions: get_box_length, get_bin_dimensions, bin_atoms_in_3D
# Req'd classes: Atoms, Neighbors
# Note: The special_neighs feature relies on atomIDs being sequential in the traj file!! 
def make_neighbor_list(r_wrap, box_bounds, RCUT, frame, special_neighs=3):
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
    
    Output: beads: Instance of Atoms class
                   beads.id2neighs[some_atomID]: 1D list of Neighbor objects containing all information about each atom's neigh list 
                                    beads.id2neighs[some_atomID].atom_ids: atom_ids in some_atomID's neigh list (1D list)
                                    beads.id2neighs[some_atomID].distances: distances for each atom ID in neigh list (1D list)
                                    beads.id2neighs[some_atomID].mol_ids: mol ids for each atom ID in neigh list (1D list)
    """
    
    number_atoms = np.shape(r_wrap)[1] - 1
    Lx, Ly, Lz = get_box_length(box_bounds, frame)
    num_xbins, num_ybins, num_zbins, _, _, _ = get_bin_dimensions(Lx, Ly, Lz, RCUT, frame)
    print "Started 3D binning procedure..."   
    num_bins = [num_xbins, num_ybins, num_zbins]
    bin2id, id2bin = bin_atoms_in_3D(num_bins, box_bounds, r_wrap, frame)
    print "Done binning!"

    print "Building neighbor list..."
    
    beads = Atoms(number_atoms)
    #id2neighbors = [ [] for n in xrange(number_atoms + 1) ]
    
    # loop over beads and their local area to build a neighbor list
    for atom_i in xrange(1, number_atoms + 1):
        neigh_list_i = Neighbors(atom_i)  # place to store neighbor info (atomID/molID/distance)
        
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
                        neigh_list_j = Neighbors(test_id_j)  # place to store neighbor info (atomID/molID/distance)
                        # if it's already in the nhood, don't bother
                        if test_id_j in neigh_list_i.atom_ids:
                            continue
                        if test_id_j == atom_i:  # Check that the two atoms aren't the same
                            continue
                        if is_on_same_molecule(atom_i, test_id_j, id2mol):  # Are atoms on the same molecule?
                            if atomIDs_too_close(atom_i, test_id_j, special_neighs):  # Are there at least special_neighs beads between the 2 atoms
                                continue

                        dx = r_wrap[frame][test_id_j][x] + xshift - r_wrap[frame][atom_i][x]
                        dy = r_wrap[frame][test_id_j][y] + yshift - r_wrap[frame][atom_i][y]
                        dz = r_wrap[frame][test_id_j][z] + zshift - r_wrap[frame][atom_i][z]
                        dr2 = dx*dx + dy*dy + dz*dz

                        # if close enough add to the neighbor list
                        if dr2 < RCUT**2:
                            neigh_list_i.add_neighbor(test_id_j, id2mol[test_id_j], np.sqrt(dr2))
                            neigh_list_j.add_neighbor(atom_i, id2mol[atom_i], np.sqrt(dr2))
                            #id2neighbors[atom_i].append(test_id_j)
                            #id2neighbors[test_id_j].append(atom_i)
        beads.store_neighs_of(atom_i, neigh_list_i)
    return beads  #return instance of Atoms class containing list of Neighbor objects w/ atomIDs, molIDs, and distances for neighbors

