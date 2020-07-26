# pppmd2
Private copy of pppmd2 - a package focused on post-processing of MD polymer simulations. Originated in Hall research group at OSU.

# Initial release of PPPMD2 (v0.0.0)
**Post-Processing Polymer Molecular Dynamics**: Suite of tools to do post-processing calculations relevant to polymers.
The primary difference between PPPMD (1) and PPPMD2 is the inclusion of the nliesen subpackage, which includes "dump_tools.py". This package adds the following functionality to the package...

**Release Date: 7/26/20**

*The original PPPMD package is included as a subpackage under ```jbrown/pppmd```*

# Initial
## nliesen/dump_tools.py
1. ```read_lammpstrj_plus```: Allows you to read in both wrapped and unwrapped coordinates from a lammps trajectory file

2. ```mini_read_lammpstrj```: Allows you to resume reading a lammps trajectory file from where you left off, without rebuilding lists such as id2type, id2mol, and mol2ids, which take a lot of time. This functionality is added with the f.seek() and f.tell() functions and enables reading in only a subset of the total frames quickly. This can be beneficial for memory usage in your program and was developed to go with a trajectory averaging script I was writing.

*Functions (1) and (2) can accept the below coordinate styles*
    
   a. wrapped (x, y, z, ix, iy, iz)  ```coordinate_type = x```
    
   b. unwrapped (xu, yu, zu)  ```coordinates_type = xu```
    
   c. scaled and wrapped (xs, ys, zs, ix, iy, iz)  ```coordinate_type = xs```

3. ```unwrap_coords```: Allows you to convert wrapped coordinates (x, y, z), image flags (ix, iy, iz), & box bounds --> unwrapped coordinates (xu, yu, zu).

4. ```wrap_coords```: Converts unwrapped coordinates (xu, yu, zu) & box bounds --> Wrapped coordinates (x, y, z) & image flags
    a. Requires function ```get_box_len```

5. ```scale_coords```: Converts wrapped coordinates (x, y, z) & box bounds --> wrapped and scaled coordinates (xs, ys, zs, ix, iy, iz)

6. ```scale_unwrapped_coords```: Converts unwrapped coords (xu, yu, zu) & box bounds --> unwrapped and scaled coordinates (xsu, ysu, zsu)

7. ```correct4_center_mass```: Resets center of mass in all frames to that of frame 0. Function accepts only unwrapped coordinates.
    *Must use with ```unwrap_coords``` and ```wrap_coords``` if you start with wrapped coordinates.*

8. ```write_lammpstrj```: Function which accepts any of the below listed forms of coordinates and "coordinate_type" and writes a file in the default style of lammps trajectory/dump files.
    a. wrapped (x, y, z, ix, iy, iz)  ```coordinate_type = x```

    b. unwrapped (xu, yu, zu)  ```coordinates_type = xu```

    c. scaled and wrapped (xs, ys, zs, ix, iy, iz)  ```coordinate_type = xs```

    d. scaled and unwrapped (xsu, ysu, zsu)  ```coordinate_type = xsu```

### Basic arrays/lists and syntax
All functions follow the procedural coding style of pppmd and work with the same kind of numpy arrays and lists present in the original code. These arrays and lists reproduced below for the sake of clarity, although detailed descriptions can be found in the code as docstrings and comments.

Position: ```r[frame, atomID, dimension]```  -- contains coordinates (3D numpy array)

Image flags: ```ir[frame, atomID, dimension]``` -- contains periodic image flags as integers  (3D numpy array) or None if unwrapped and ```flags_when_unwrap=False```

Box dimensions: ```box_bounds[frame, dimension, low(0)/high(1)]``` -- contains low/high x, y, and z bounds  (3D numpy array)

Atom types: ```id2type[atomID]``` -- corresponding integer-valued atom type for each atomID  (1D numpy array)

AtomID-->molID: ```id2mol[atomID]``` -- molecule/moleculeID for chain atom/atomID belongs to  (1D numpy array)

Frame-->Timestep: ```timestep[frame]``` -- contains timestep of frame (1D numpy array)

MolD-->atomIDs```mol2ids``` -- List whose entries are numpy arrays containing all atomIDs for a selected molID. The list is indexed by molID. (list of 1D np arrays)

#### jbrown/pppmd

1. ```read_lammpstrj```: Allows you to read in only wrapped coordinates (scaled or unscaled) from a lammps trajectory file.

2. ```MSD```: Calculates the mean square displacement for each bead type using the wrapped coordinates, image flags, box bounds, and atom types. Output is a dictionary of numpy arrays, where the numpy arrays are indexed by frame.

    *does NOT do any sort of block averaging*

    *assumes mass = 1 for all beads*

    *does not account for changes in box size*

3. ```gofr```: Computes radial distribution function using the unscaled, but wrapped coordinates (```r```), the box bounds (```box_bounds```), and chosen
bin size (```bin_size```) in unscaled units of distance.
     *does not account for changes in box size*
     *does not distinguish between atoms on the same molecule and other atoms*

4. ```Sofk```: Uses unscaled coordinates, box bounds, and number of bins to directly calculate sofk and the partial structure factors between all pairs of atom types, rather than using a fourier transform of h(r)=g(r)-1. 

     *does not account for changes in box size*

     *does not distinguish between atoms on the same molecule and other atoms*

5. ```end2end_autocorr```: Converts unscaled (but wrapped) coords, image flags, box bounds, and the list ```mol2ids```, which contains the atomIDs belonging to each molecule (the list is indexed by molID), to an end to end autocorrelation (ACF) function. The outputted ACF is a 1D array indexed by frame count. 

    *all the listings in mol2ids will be used and averaged together*

    *it is assumed that the end-to-end vector is the one between the lowest and highest id in each molecule (if this is not the case, you'd have to mess with mol2ids, e.g. make it only contain the ids of the two end beads)*

    *scaled by the average end-to-end vector at frame 0, so that e2e_autocorr[0]=1.0*

##### Note on semantics for versioning:
Much of our code will rely on these scripts. To avoid breaking code (i.e. have backwards compatible code) I will follow the scheme outlined below when adding "releases" to keep track of our code as it changes. 

**Version schematics: MAJOR.MINOR.PATCH**

    **MAJOR version**: Functions no longer return the expected objects/variables in the ways that our code relied on

         *e.g. ```read_lammpstrj_plus``` no longer accepts returns objects in the expected order which breaks backwards compatibility*

    **MINOR version**: when you add functionality in a backwards compatible manner, and

        *e.g. when you add a function or a subpackage, without affecting how other functions accept or return objects, or when extend the functionality without breaking how it accepts/returns objects.*

    **PATCH version**: when you make backwards compatible bug fixes.

        *e.g. fixing a bug in the code without affecting how functions accept/return things.

For more details see [this page](https://semver.org/).
