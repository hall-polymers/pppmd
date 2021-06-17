#%%

import numpy as np
import sys
from math import floor
sys.path.append("../")
import dump_tools as dtools

#%%

def find_beads_of_type(bead_type, id2type):
    """Takes in the type of an atom, and the array which stores atom types (indexed by atomID)
    and generates a list of atomIDs corresponding to the selected type."""
    print "Identifying beads of type "+str(bead_type)+" and outputting relevant atomIDs."
    atomIDs = np.where(id2type==bead_type)[0]
    return atomIDs


def net_mass_beads(type_bead, mass, id_to_type):
    total_mass = 0
    print "Finding beads of type "+str(type_bead)
    atomIDs = find_beads_of_type(type_bead, id_to_type)
    number_beads = len(atomIDs)
    print "number of beads of this type is "+str(number_beads)
    mass_beads = number_beads*mass
    print "net mass of these beads is "+str(mass_beads)

    return mass_beads, atomIDs


# Add ability to subtract out center of mass drift from coordinates
# NOTE: By default assume x, y, and z should all be corrected
# NOTE: The below comes from the PGN bond vector field scripts
def correct4_center_mass_plus(r_unwrap, id_2_type, type_2_mass, x_is_periodic = True, y_is_periodic = True, z_is_periodic = True):
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
            mr_t[atomIDs_of_type, :] = r_t[atomIDs_of_type, :]*mass

        com_position_t = np.sum(mr_t, axis = 0)/total_mass  # Center-of-mass position for frame t
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
        #if x_is_periodic and y_is_periodic and z_is_periodic:
        #    r_corrected[t, :, :] = r_unwrap[t, :, :] - change_com  # works via broadcasting - subtract change in COM
        if x_is_periodic:
            r_corrected[t, :, x] = r_unwrap[t, :, x] - change_com[x]  # works via broadcasting - subtract change in X component of COM
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


# %%

# Written by NLiesen during our Summer @ AFRL/RX: 6/9/21
# Provides the ability to average or smooth trajectories a la VMD
def smooth_lammpstrj(input_lammpstrj, nFrames, windowSize, windowOffset, coords_are_unwrapped, skip_frame_zero = False,
 is_com_correction_required = False, wrap_before_writing = False, type2mass = None, is_periodic = (True, True, True),
 output_fname = 'smoothed.lammpstrj'):
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
    
    if type(skip_frame_zero) != bool or type(is_com_correction_required) != bool\
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
        nWindows = int(floor( ((nFrames - 1.) - windowSize) / windowOffset) + 1) 
    else:
        nWindows = int(floor( float(nFrames - windowSize) / windowOffset) + 1)
    
    print "Total number of windows: {}".format(nWindows)

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

        print "Skipping {0} frames for window {1}".format(skip_nFrames, windowID)

        if windowID == 0:  # read_lammpstrj_plus is smart enough to detect whether coords are wrapped
            rw, _, tsteps, boxbds_window, id2type, id2mol, _ = dtools.read_lammpstrj_plus(fname = input_lammpstrj,
             num_frames = windowSize, skip_beginning = skip_nFrames, skip_between = 0, flags_when_unwrap = False)
        elif windowID > 0:
            rw, _, tsteps, boxbds_window, _, _, _ = dtools.read_lammpstrj_plus(fname = input_lammpstrj,
             num_frames = windowSize, skip_beginning = skip_nFrames, skip_between = 0, flags_when_unwrap = False)

        print "Finished reading coordinates for window {}".format(windowID)

        #windowTime = (np.max(tsteps) - np.min(tsteps)) / 2.

        # average these coordinates to get the average positions of the atoms within the window
        if coords_are_unwrapped:
            if is_com_correction_required:
                # first we need to remove the center of mass motion
                rw = correct4_center_mass_plus(rw, id2type, type2mass, x_is_periodic = is_periodic[x],
                 y_is_periodic = is_periodic[y], z_is_periodic = is_periodic[z])

            rw_ave = np.mean(rw, axis = timeAxis)  # average atomic position vectors over time
            rw_ave = rw_ave.reshape((1, rw.shape[atomAxis], rw.shape[dimAxis]))
            boxbds_ave = np.mean(boxbds_window, axis = timeAxis)  # average lower and upper bounds of box over time
            boxbds_ave = boxbds_ave.reshape((1, boxbds_window.shape[1], boxbds_window.shape[2]))
            if wrap_before_writing:
                # Now we will want to re-wrap these coordinates
                rw_ave, ir_ave = dtools.wrap_coords(rw_ave, boxbds_ave)
                # Now that we have processed the required coordinates for each window we must write to a .lammpstrj
                # file.
                dtools.write_lammpstrj(output_fname, boxbds_ave,
                 np.array([windowID]), id2mol, id2type, rw_ave, image_flags = ir_ave,
                 coordinate_type = 'x', write_mode = writeMode )
            else:
                dtools.write_lammpstrj(output_fname, boxbds_ave,
                 np.array([windowID]), id2mol, id2type, rw_ave, coordinate_type = 'xu', write_mode = writeMode )
        else:  # TODO: Add ability to read in wrapped coordinates using pppmd utilities
            raise Exception('Coordinates must be unwrapped prior to smoothing')

        print "Finished window number: {}".format(windowID)

#%%

if __name__ == '__main__':
    # non-overlapping windows for test
    number_frames = 101
    window_size = 100 # Each window consists of 10 frames
    window_offset = 100

    myfile = 'example_file.lammpstrj'

    type_to_mass = {  # dictionary containing atom type: corresponding mass
        1: 1.,
        2: 1.,
        3: 445.1
    }

    smooth_lammpstrj(myfile, number_frames, window_size, window_offset, coords_are_unwrapped = True, skip_frame_zero = True,
     is_com_correction_required = True, wrap_before_writing = False, type2mass = type_to_mass, is_periodic = (True, True, False),
     output_fname = 'smoothed_over_100frames.lammpstrj')

# %%
