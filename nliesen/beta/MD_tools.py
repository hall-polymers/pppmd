#!/usr/bin/env python2
# coding: utf-8

#%%

# Purpose: Original focus on getting neighbors within some spherical cutoff of each atom
# Author: N. Liesen(NTL)
# Resources: Most of the code is stolen from cluster_v3.1.py, which was a collaborative
# script from the Hall group (Jon/Anne/Kyaw/Janani/Conner and others?) and included the
# main code required for neighbor listing
# Date Created: 6/27/20 -- NTL
# Last Modified: 8/20/20 -- NTL

import numpy as np
import random
from math import pi, sqrt, cos, sin, floor

import sys
sys.path.append("../")
import dump_tools as dtools

import matplotlib.style
import matplotlib.font_manager
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
plt.interactive(True)
from mpl_toolkits.mplot3d import Axes3D

import pprint as pp

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
    vx = sqrt(1-vz**2)*cos(theta)
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

class Contacts:
    """
    This class will focus in on close contacts and their positions
    """
    def __init__(self, current_frame, contact_cutoff):
        self.contact_pairs = []  # list of (atomi, atomj) sets for contacts; sets used b/c order
        # of atoms doesn't matter
        self.chain_pairs = []
        self.contact_positions = []  # list of position vectors (arrays) for contacts
        self.number_contacts = 0
        self.frame = current_frame
        self.contact_rcut = contact_cutoff
        self.is_wrapped = False
        self.contact_flags = None
        self.contactID2cluster = {}  # list of integers identifying the cluster to which each contact
        # belongs. The array is indexed by contactID, an integer uniquely identifying each contact.
    
    def add_contact(self, ID_atomI, ID_atomJ, r_atomI, r_atomJ, ID_molA, ID_molB, check_for_duplicate = True):
        """
        Checks if the i/j pair of atoms already has been counted - if not, stores the pair of atoms
        as a contact along with the average position vector of the atoms making up the contact.
        """

        if (not check_for_duplicate) or ({ID_atomI, ID_atomJ} not in self.contact_pairs):  # or is short circuit
            self.contact_pairs.append({ID_atomI, ID_atomJ})  # add {i, j} pair to list
            self.chain_pairs.append({ID_molA, ID_molB})
            self.contact_positions.append((r_atomI + r_atomJ) / 2.)  # add ave contact position vector
            self.number_contacts += 1

    def wrap_contacts(self, boxbounds):
        """
        The way we calculate average positions for contacts, they must be re-wrapped; this is b/c some
        contact positions will fall between the periodic image of one chain, and the main image of another.
        Wrapping will ensure that all identified contact's position vectors will fall within the central
        simulation box.
        """
        r_contacts = np.array(self.contact_positions).reshape(1, len(self.contact_positions), 3)
        self.contact_positions, self.contact_flags = dtools.wrap_coords(r_contacts, boxbounds)
        self.is_wrapped = True

    def contacts_by_chain(self):
        """
        NOTE: Cannot use normal sets as keys, dict keys must be immutable
        In the end I decided to use frozen sets as the keys (these are just
        normal sets, but immutable). I wanted to use sets b/c the order in which
        the ij chain pairs are passed as keys shouldn't matter. Keys ij and ji
        should result in the same corresponding list of contactIDs.
        """
        chains2contacts = {}
        for contactID in range(self.number_contacts):
            chain_i, chain_j = self.chain_pairs[contactID]
            chain_ij = frozenset((chain_i, chain_j))
            if chain_ij in chains2contacts.keys():
                chains2contacts[chain_ij].append(contactID)
            else:
                chains2contacts[chain_ij] = [contactID]
        self.chains_to_contacts = chains2contacts 

    def get_nBeads_separating_contacts(self, contactIDx, contactIDy, atomID2mol):
        """
        This function will accept two contactIDs, provided that the two contacts belong to the
        same ij pair of chains (i.e. contact IDs/points x and y both occur as a result of contact
        between beads on chains i and j); the function will identify the number of beads
        separating bead B_ix (the bead on chain i participating in contact point x) and B_iy
        (the bead on chain participating in contact point y). It will similarly find the number
        of beads separating beads B_jx and B_jy. The returned value will correspond to
        min((B_ix - B_iy - 1), (B_jx - B_jy - 1)). i.e. if the number of beads separating bead
        B_ix and B_iy is 5, and the number of beads separating beads B_jx and B_jy is 15 then
        we will return 5; the justification being that if any bead in contact x is w/in some number
        of beads (usually 20) of a bead participating in contact y, then the two contacts belong to
        the same cluster. 

        The purpose of returning the number of monomers separating these pairs of beads on chains
        i and j is to determine whether contacts x and y belong to the same cluster or not. The
        clustering will likely need to be done recursively...
        """

        # First, check that the two contact points actually originate from the same two chains
        if not self.chain_pairs[contactIDx] == self.chain_pairs[contactIDy]:
            raise Exception('The 2 contact points specified originate from different pairs of chains')
        elif contactIDx == contactIDy:
            raise Exception('The two contact IDs are the same.')
        # NOTE: Assume that a lower integer-valued moleculeID implies a lower integer-valued atomID

        atomIDs = {}
        chain_i, chain_j = sorted(list(self.chain_pairs[contactIDx]))
        for contact in (contactIDx, contactIDy):
            for atom in self.contact_pairs[contact]:
                if id2mol[atom] == chain_i:
                    atomIDs[(contact, chain_i)] = atom
                elif id2mol[atom] == chain_j:
                    atomIDs[(contact, chain_j)] = atom

        beads_btwn_contacts_i = abs(atomIDs[(contactIDx, chain_i)] - atomIDs[(contactIDy, chain_i)])
        beads_btwn_contacts_j = abs(atomIDs[(contactIDx, chain_j)] - atomIDs[(contactIDy, chain_j)])

        return min(beads_btwn_contacts_i, beads_btwn_contacts_j)

    def __str__(self):
        return "In frame {0} there are {1} close contacts, given a distance of {2}".format(self.frame,
        self.number_contacts, self.contact_rcut)

class Neighbors:
    """Each neighbor will have it's atomID stored in self.atom_ids. The atomID's corresponding molID will be stored in
    self.mol_ids. The distance used to determine that the atom belongs in this neighbor list is stored in self.distances.
    The add_neighbor method will be called to add the atomID, molID, and distance into the proper list attribute."""
    def __init__(self,  central_atomID):
        self.center_atomID = central_atomID  # To verify the neighbor list is 4 correct atom
        self.atom_ids = []
        self.mol_ids = []
        self.distances = []

        self.nneigh = None  # atomID of nearest neighbor
        self.nneigh_dist = None  # distance central atom --> nearest neighbor
        self.farthest_neigh = None  # atomID of farthest neighbor
        self.farthest_dist = None  # distance central atom --> farthest neighbor

    def add_neighbor(self, new_atomID, new_molID, dist):
        self.atom_ids.append(new_atomID)
        self.mol_ids.append(new_molID)
        self.distances.append(dist)
        return

    def get_nearest_neighbor(self):
        """
        This method will find the neighbor nearest to the central atom.

        returns:
        nneigh: the atomID of the nearest neighbor
        nneigh_dist: the distance between the central atom and this nearest neighbor
        """
        number_neighbors = len(self.atom_ids)
        if number_neighbors > 0:
            nneigh = self.atom_ids[np.argmin(np.array(self.distances))]
            nneigh_dist = np.min(np.array(self.distances))
            self.nneigh = nneigh
            self.nneigh_dist = nneigh_dist
            return nneigh, nneigh_dist
        elif number_neighbors == 0:
            print("no neighbors!")
            return None, None

    def get_farthest_neighbor(self):
        """
        This method will find the neighbor farthest from the central atom.

        returns:
        farthest_neigh: the atomID of the neighbor farthest from the central atom
        farthest_dist: the distance between the central atom and the farthest neighbor
        """
        number_neighbors = len(self.atom_ids)
        if number_neighbors > 0:
            farthest_neigh = self.atom_ids[np.argmin(np.array(self.distances))]
            farthest_dist = np.max(np.array(self.distances))
            self.farthest_neigh = farthest_neigh
            self.farthest_dist = farthest_dist
            return farthest_neigh, farthest_dist
        elif number_neighbors == 0:
            print("no neighbors!")
            return None, None

    def __str__(self):  # modify print behavior of class instances
        number_neighbors = len(self.atom_ids)
        if number_neighbors > 0:
            #nneigh_dist = np.min(np.array(self.distances))
            #farthest_dist = np.max(np.array(self.distances))
            _, nneighbor_dist = self.get_nearest_neighbor()
            _, farthest_dist = self.get_farthest_neighbor()
            return "atomID {0} has {1} neighbors. ".format(self.center_atomID, number_neighbors) + \
                "The nearest and farthest neighbors are at a distance of {0} and {1} respectively.".format(nneighbor_dist, farthest_dist)
        elif number_neighbors == 0:
            return "atomid {0} has no neighbors".format(self.center_atomID)

#%%

# Req'd functions: get_box_length, get_bin_dimensions, bin_atoms_in_3D
# Req'd classes: Neighbors
# Note: The special_neighs feature relies on atomIDs being sequential in the traj file!! 
def make_neighbor_list(r_wrap, box_bounds, RCUT, frame, atom2molid, id2type, special_neighs=3,
 excluded_types = [], is_periodic = (True, True, True)):
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


    # argument handling
    if (type(r_wrap) != np.ndarray) or (r_wrap.ndim != 3):
        raise TypeError('Coordinates, r[frame, atomID, dim], must be passed as a 3D numpy array')
    elif (type(box_bounds) != np.ndarray) or (box_bounds.ndim != 3):
        raise TypeError('Box Bounds, box_bounds[frame, dim, lower/upper] must be passed as 3D numpy array')

    number_atoms = np.shape(r_wrap)[1] - 1  # req'd for handling of the remaining arguments

    if type(frame) != int:
        raise TypeError('frame must be integer-valued')
    elif type(special_neighs) != int:
        raise TypeError('special neighbors parameter must be integer-valued')

    if type(excluded_types) != list:
        raise TypeError('The excluded bead types must be passed as a list')

    #if type(id2type) != np.ndarray or (len(id2type) != (number_atoms + 1))\
    # or (id2type.ndim != 1):
    #    raise TypeError('The atomID --> atomtype array, id2type, must be a 1D numpy array whose length is\
    #    commensurate with the number of atoms')
    #elif type(atom2molid) != np.ndarray or (len(atom2molid) != (number_atoms + 1))\
    # or (atom2molid.ndim != 1):
    #    raise TypeError('The atomID --> molid array, atom2molid, must be a 1D numpy array whose length is\
    #    commensurate with the number of atoms')

    if type(is_periodic) != tuple or len(is_periodic) != 3:
        raise TypeError('is_periodic must be passed as a tuple of length 3')

    for xyz_flag in is_periodic:
        if type(xyz_flag) != bool:
            raise TypeError('Boundary condition flags in is_periodic must be boolean')
    
    # aliases
    x = 0
    y = 1
    z = 2

    Lx, Ly, Lz = get_box_length(box_bounds, frame)
    num_xbins, num_ybins, num_zbins, _, _, _ = get_bin_dimensions(Lx, Ly, Lz, RCUT, frame)
    print "Started 3D binning procedure..."
    num_bins = [num_xbins, num_ybins, num_zbins]
    bin2id, id2bin = bin_atoms_in_3D(num_bins, box_bounds, r_wrap, frame)
    print "Done binning!"

    print "Building neighbor list..."

    id2neighbors = [ Neighbors(atom_i) for atom_i in range(0, number_atoms + 1)]  # List of neighbor objects
    id2neighbors[0] = None  # atomIDs are 1-indexed
    my_contacts = Contacts(frame, RCUT)

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
        """
        We must remember that our PGN systems aren't periodic along the z-dimension. This means
        that near the edges of the monolayer (meaning the +z or -z boundary) we need to check
        fewer bins, since there is just a wall or vacuum aloong the -z and +z boundaries respectively.
        To handle cases where the boundary conditions are fixed/shrink-wrapped/or some other
        non-periodic situation I've added...
        i. An argument indicating whether each dimension is periodic or not.
        ii. A check that the relevant dimension is periodic when deciding to check a periodic image of a
        voxel for neighbors.
        """
        for i in xrange(atoms_xbin-1, atoms_xbin+2):
            # if not periodic
            if i == -1 and not is_periodic[x]:
                continue
            elif i == last_xbin + 1 and not is_periodic[x]:
                continue
            # if periodic
            elif i == -1 and is_periodic[x]:
                i = last_xbin  # Periodic
                xshift = -Lx
            elif i == last_xbin + 1 and is_periodic[x]:
                i = first_bin  # Periodic
                xshift = Lx
            else:
                xshift = 0

            for j in xrange(atoms_ybin-1, atoms_ybin+2):
                # if not periodic
                if j == -1 and not is_periodic[y]:
                    continue 
                elif j == last_ybin + 1 and not is_periodic[y]:
                    continue
                # if periodic
                elif j == -1 and is_periodic[y]:
                    j = last_ybin
                    yshift = -Ly
                elif j == last_ybin + 1 and is_periodic[y]:
                    j = first_bin
                    yshift = Ly
                else:
                    yshift = 0

                for k in xrange(atoms_zbin-1, atoms_zbin+2):
                    # if not periodic
                    if k == -1 and not is_periodic[z]:
                        continue
                    elif k == last_zbin + 1 and not is_periodic[z]:
                        continue
                    # if periodic
                    elif k == -1 and is_periodic[z]:
                        k = last_zbin
                        zshift = -Lz
                    elif k == last_zbin + 1 and is_periodic[z]:
                        k = first_bin
                        zshift = Lz
                    else:
                        zshift = 0

                    # loop over the beads in this box and calculate the distance
                    for test_id_j in bin2id[i][j][k]:
                        if id2type[atom_i] in excluded_types:
                            continue 
                        if id2type[test_id_j] in excluded_types:
                            continue
                        # if it's already in the nhood, don't bother
                        if test_id_j in id2neighbors[atom_i].atom_ids:  # check if test atom is already in atom i's nlist
                            continue
                        if test_id_j == atom_i:  # Check that the two atoms aren't the same
                            continue
                        if is_on_same_molecule(atom_i, test_id_j, atom2molid):  # Are atoms on the same molecule?
                            if atomIDs_too_close(atom_i, test_id_j, special_neighs):  # Are there at least special_neighs beads between the 2 atoms
                                continue

                        # These contact positions will need to be re-wrapped
                        rx_atom_j = r_wrap[frame][test_id_j][x] + xshift
                        rx_atom_i = r_wrap[frame][atom_i][x]
                        ry_atom_j = r_wrap[frame][test_id_j][y] + yshift
                        ry_atom_i = r_wrap[frame][atom_i][y]
                        rz_atom_j = r_wrap[frame][test_id_j][z] + zshift
                        rz_atom_i = r_wrap[frame][atom_i][z]

                        #dx = r_wrap[frame][test_id_j][x] + xshift - r_wrap[frame][atom_i][x]
                        #dy = r_wrap[frame][test_id_j][y] + yshift - r_wrap[frame][atom_i][y]
                        #dz = r_wrap[frame][test_id_j][z] + zshift - r_wrap[frame][atom_i][z]
                        dx = rx_atom_j - rx_atom_i
                        dy = ry_atom_j - ry_atom_i
                        dz = rz_atom_j - rz_atom_i
                        dr2 = dx*dx + dy*dy + dz*dz

                        # if close enough add to the neighbor list
                        if dr2 < RCUT**2:
                            id2neighbors[atom_i].add_neighbor(test_id_j, atom2molid[test_id_j], np.sqrt(dr2))
                            id2neighbors[test_id_j].add_neighbor(atom_i, atom2molid[atom_i], np.sqrt(dr2))
                            r_atom_i = np.array([rx_atom_i, ry_atom_i, rz_atom_i])
                            r_atom_j = np.array([rx_atom_j, ry_atom_j, rz_atom_j])
                            my_contacts.add_contact(atom_i, test_id_j, r_atom_i, r_atom_j,
                            atom2molid[atom_i], atom2molid[test_id_j], check_for_duplicate = False)  # store as close contact
    my_contacts.wrap_contacts(box_bounds)  # some contacts will end up outside the box -- easiest way to fix is to just re-wrap
    return id2neighbors, my_contacts  #return instance of Atoms class containing list of Neighbor objects w/ atomIDs, molIDs, and distances for neighbors


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
    id2type = [0]

    # generate atomic positions such that the atoms fall onto the surface
    # of spheres of different radii
    frame = 0
    atomid = 1 
    r[frame, atomid] = my_center
    id2mol.append(atomid)
    id2type.append(atomid)
    atomid += 1
    for radius in my_radii:
        for ii in range(points_per_radius):
            my_vector = make_random_vector(desired_radius = radius)
            r[frame, atomid] = my_center + my_vector 
            id2mol.append(atomid)  # new molecule for every placed bead
            id2type.append(atomid)  # new type for every placed bead
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
    id2neighlist, _ = make_neighbor_list(r, boxbds, 4, frame, id2mol, id2type, special_neighs=0)

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

#%%

# Now running tests with real lammpstrj files

# read in smoothed lammpstrj file
mylammpstrj = 'smoothed_over_10frames.lammpstrj'  # these coordinates are unwrapped
r, _, tsteps, boxbds, id2type, id2mol, mol2ids = dtools.read_lammpstrj_plus(fname = mylammpstrj,
 flags_when_unwrap = False)

# wrap smoothed coordinates
r, ir = wrap_coords(r, boxbds)

chain_length = 160
number_frames = r.shape[0]

my_write_mode = 'w'
#rCut = 0.5  # sigma
#for frame in range(number_frames):
#rCut_list = np.arange(0.2, 1.3, 0.1)
rCut_list = [0.3, 0.4]

nContacts = []
for rCutoff in rCut_list:
    for frame in range(1):
        # get neighbors using smoothed and wrapped coordinates
        id2neighlist, contacts = make_neighbor_list(r, boxbds, rCutoff, frame, id2mol, id2type, special_neighs = chain_length,
         excluded_types = [3], is_periodic = (True, True, False))  # z isn't periodic

        # writing our contacts to a file...
        contactID2mol = [1 for i in range(contacts.number_contacts)]
        contactID2type = [1 for i in range(contacts.number_contacts)]
        my_box_bounds = boxbds[frame].reshape(1, 3, 2)

        dtools.write_lammpstrj('contacts_' + str(rCutoff) + 'sigma_first_frame.lammpstrj', my_box_bounds,
         np.array([frame]), contactID2mol, contactID2type, contacts.contact_positions, image_flags = contacts.contact_flags,
         coordinate_type = 'x', write_mode = my_write_mode) 

        #dtools.write_lammpstrj('wrapped_smoothed_first_frame.lammpstrj', my_box_bounds,
        # np.array([frame]), id2mol, id2type, r[frame, :, :].reshape(1, r.shape[1], r.shape[2]),
        #  image_flags = ir[frame, :, :].reshape(1, ir.shape[1], ir.shape[2]),
        # coordinate_type = 'x', write_mode = my_write_mode) 

        nContacts.append(contacts.number_contacts)

        my_write_mode = 'a'  # append to lammpstrj files we just created for the first frame

# Now that we have managed to generate lammpstrj files containing all the contacts we should
# generate a heat map using the positions of these atoms. One option is to read the .lammpstrj
# file and then turn the results into a heat map. We can use the code for binning bond vectors
# from our voxelization scripts to do this in order to avoid re-writing fairly similar code.
# The basic steps involved should be to convert to scaled coordinates, re-scale based on the
# average box dimensions, and then bin in 3D using the average box dimensions. The heat map can
# be visualized using paraview. The results can be output in the same manner as the bond counts.

#%%

"""
# pick an arbitrary atomID and find the nearest contact in its neighbor list
my_atomID = 30000
print(id2neighlist[my_atomID])
neighbors = id2neighlist[my_atomID].atom_ids
neighbors.sort() 
print neighbors

for atom, distance in zip(id2neighlist[my_atomID].atom_ids, id2neighlist[my_atomID].distances):
    print "atom {0} is at a distance of {1} from atom {2}".format(atom, distance, 
    id2neighlist[my_atomID].center_atomID)

id2neighlist[my_atomID].get_nearest_neighbor()
nearest_neighbor = id2neighlist[my_atomID].nneigh
dist2nearest_neighbor = id2neighlist[my_atomID].nneigh_dist

# now find the moleculeID for this nearest neighbor
central_mol = id2mol[my_atomID]
nearest_mol = id2mol[nearest_neighbor]

print(central_mol)
print(nearest_mol)
"""

# %%

# %%

# for two chains create a list of atom pairs and the distances between the two atoms

# lets start with an arbitrary chain
chainID = 450
atoms_in_chain = mol2ids[chainID]

# for each atom in the selected chain find all atoms which are the nearest neighbor to some
# atom on the chain
atom_i = []
chain_i = []
atom_j = []
chain_j = []
neigh_distances = []
for atomID in atoms_in_chain:
    nearest_neighbor, nearest_neigh_distance = id2neighlist[atomID].get_nearest_neighbor()
    if nearest_neighbor is not None:
        atom_i.append(atomID)
        chain_i.append(id2mol[atomID])
        atom_j.append(nearest_neighbor)
        chain_j.append(id2mol[nearest_neighbor])
        neigh_distances.append(nearest_neigh_distance)

# find the chain which passes nearest to an atom on chain i
nearest_chain = chain_j[np.argmin(np.array(neigh_distances))]
print chainID
print nearest_chain

#%%

cutoffs = np.arange(0.2, 1.3, 0.1)
coeff = 316
ypts = coeff*cutoffs**3

# %%
