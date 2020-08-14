## PPPMD Ion Dynamics Package
**Developed by Kuan-Hsuan Kevin Shen**

This sub-package includes python analysis scripts for ion-containing systems, including ion pairing/clustering and ion conductivity analyses.

1. ```special_read```: Allows you to read in coordinates of certain types of beads from a lammps trajectory file (slightly modified version of original ```pppmd.read_lammpstrj``` function. 

**Input: fname, types, num_frames, skip_beginning, skip_between**
- fname: filename string or 'stdin' (or a value that evaluates to false) for reading from standard in 
- types: *list* of atom types to extract coordinates for
- num_frames: optional number of frames to read before stopping, defaults to reading in all frames
- skip_beginning: skip this many frames at the beginning of the dump file
- skip_between: skip this many frames between saved frames

**Output: r, ir, timestep, box_bounds, id2type, id2index**
- r: num_frames by num_atoms+1 by 3 array of wrapped and unscaled coordinates (indexed by frame number then atom id)
- ir: num_frames by num_atoms+1 by 3 array of image flags
- timestep: num_frames length array of timesteps
- box_bounds: 3D array to store boundaries of the box, indexed by frame, x/y/z, then lower/upper
- id2type, id2index: num_atoms+1 length arrays to map atom id to type and index id (if available, may be None)

This function also outputs the amounts of total beads and selected beads as follows
```
===============Summary===============
Total number of beads =  num_atoms
Total number of selected beads = num_types
```

2. ```buildnlist```: Builds a neighbor list for each bead and return a id2neighbors: num_frames by num_atoms+1 by 3 array of neighbor list. Optional: ```nearest=True```: only have the nearest bead in the neighbor list.


**General block averaging syntax used in these functions:**
- ```nBlock```: Number of blocks
- ```blockSize```: Number of time frames in each block
- ```intrvl```: Number of time frames between the starting point of each block

The following schematic is an example of ```nBlock=3``` that uses all the time frames of the entire simulation (in this case, ```(nBlock-1)*intrvl + blockSize = len(r)```:

```
!--------------------------Entire simulation--------------------------!
|________________blockSize________________|
    intrvl    |________________blockSize________________|
                  intrvl    |________________blockSize________________|

```

*Detailed explanation of usage of kshen/ion-dynamics is coming.*
*Example usages of kshen/ion-dynamics functions are available in [kshen/ion-dynamics/example](https://github.com/hall-polymers/pppmd2/tree/development/kshen).*
