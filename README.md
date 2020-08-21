# PPPMD development
Development copy of pppmd - a package focused on post-processing of MD polymer simulations. Originated in Hall research group at OSU.
---
## Table of Contents
- [Basic Arrays/Lists and Syntax](#basic)
- [Usage of Functions](#usage)
    - [jbrown/pppmd](https://github.com/hall-polymers/pppmd2/tree/development/jbrown)
    - [nliesen/dump_tools](https://github.com/hall-polymers/pppmd2/tree/development/nliesen)
    - [nliesen/beta/MD_tools](https://github.com/hall-polymers/pppmd2/tree/development/nliesen/beta)
    - [kshen/ion_dynamics](https://github.com/hall-polymers/pppmd2/tree/development/kshen)
- [Release Descriptions](#release) 
---
## Basic Arrays/Lists and Syntax <a name="basic"/>
All functions follow the procedural coding style of pppmd and work with the same kind of numpy arrays and lists present in the original code. These arrays and lists reproduced below for the sake of clarity, although detailed descriptions can be found in the code as docstrings and comments.

- Position: ```r[frame, atomID, dimension]```  -- contains coordinates (3D numpy array)
- Image flags: ```ir[frame, atomID, dimension]``` -- contains periodic image flags as integers  (3D numpy array) or None if unwrapped and ```flags_when_unwrap=False```
-  Dimensions: ```box_bounds[frame, dimension, low(0)/high(1)]``` -- contains low/high x, y, and z bounds  (3D numpy array)
- Atom types: ```id2type[atomID]``` -- corresponding integer-valued atom type for each atomID  (1D numpy array)
- AtomID &rarr; molID: ```id2mol[atomID]``` -- molecule/moleculeID for chain atom/atomID belongs to  (1D numpy array)
- Frame &rarr; Timestep: ```timestep[frame]``` -- contains timestep of frame (1D numpy array)
- MolD &rarr; atomIDs```mol2ids``` -- List whose entries are numpy arrays containing all atomIDs for a selected molID. The list is indexed by molID. (list of 1D np arrays)


## Usage of Functions <a name="usage"/>
We describe all the details of the functions in the subfolders:
- [jbrown/pppmd](https://github.com/hall-polymers/pppmd/tree/development/jbrown) -- PPPMD package that includes basic functions to read in LAMMPS trajectory files, calculate particle mean square displacement, radial distribution function, structure factor, and end-to-end autocorrelation function.
- [nliesen/dump_tools](https://github.com/hall-polymers/pppmd/tree/development/nliesen) -- Advanced PPPMD package that allows 1) read in both wrapped and unwrapped coordinates from a lammps trajectory file and 2) resume reading a lammps trajectory file from where you left off.
- [nliesen/beta/MD_tools](https://github.com/hall-polymers/pppmd/tree/development/nliesen/beta) -- Recently added package with the ability to bin atoms into 3D cells as their trajectories evolve over time, and to generate neighbor lists for each atom using this binning procedure to speed up the neighbor listing by only checking nearby cells.
- [kshen/ion_dynamics](https://github.com/hall-polymers/pppmd/tree/development/kshen) -- PPPMD Ion Dynamics Package that includes functions to analyze ion pairing/clustering and calculate ion conductivity.

## Release Descriptions <a name="release"/>

### Initial release of PPPMD (v0.1.0)
**Post-Processing Polymer Molecular Dynamics**: Suite of tools to do post-processing calculations relevant to polymers.
The primary changes in made in this release is the inclusion of the nliesen subpackage, which includes "dump_tools.py", which is described in more detail [here](https://github.com/hall-polymers/pppmd/tree/development/nliesen).

**Release Date: 7/26/20**

*The original functions that made up the PPPMD package (prior to being posted on github) is included as a subpackage under ```jbrown/pppmd```*


##### Note on semantics for versioning:
Much of our code will rely on these scripts. To avoid breaking code (i.e. have backwards compatible code) I will follow the scheme outlined below when adding "releases" to keep track of our code as it changes. 

**Version semantics: MAJOR.MINOR.PATCH**

    **MAJOR version**: Functions no longer return the expected objects/variables in the ways that our code relied on

         *e.g. ```read_lammpstrj_plus``` no longer accepts returns objects in the expected order which breaks backwards compatibility*

    **MINOR version**: when you add functionality in a backwards compatible manner, and

        *e.g. when you add a function or a subpackage, without affecting how other functions accept or return objects, or when extend the functionality without breaking how it accepts/returns objects.*

    **PATCH version**: when you make backwards compatible bug fixes.

        *e.g. fixing a bug in the code without affecting how functions accept/return things.

For more details see [this page](https://semver.org/).
