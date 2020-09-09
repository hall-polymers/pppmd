## PPPMD Ion Dynamics Package
**Developed by Kuan-Hsuan Kevin Shen**

This sub-package includes python analysis scripts for ion-containing systems, including ion pairing/clustering and ion conductivity analyses.

| Methods                                                                                   | Description                                                                                                            |
|-------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------|
| [special_read](#special_read)(fname, types[, num_frames, skip_beginning, skip_between])   | Read in coordinates of certain types of beads from a LAMMPS trajectory file                                            |
| [wrap](#wrap)(r, box_bounds)                                                              | Simply wrap coordinates                                                                                                |
| [binbox](#binbox)(r, boxsize, bound_lo, dist_range)                                       | Bin the box and assigns each atom into the corresponding bin                                                          |
| [buildnlist](#buildnlist)(r, bin2id, id2bin, bins, boxsize, id2type, dist_range, nearest) | Build a neighbor list for each bead and return a id2neighbors |
| [ipcorr](#ipcorr)(r, id2neighbors) | Calculate ion pair autocorrelation function|
| [cluster_autocorr](#cluster_autocorr)(frame, sflag) | Calculate ion cluster autocorrelation function|
| [buildclusterlist](#buildclusterlist)(r, id2neighbors_1stshell, id2neighbors_2ndshell) | Calculate average ion cluster size and ion cluster autocorrelation function, and analyze local ion morphology (categorize ions into 4 groups: aggregate, contact ion pair, monomer-separated ion pair, free ion)|

1. ```special_read```: Allows you to read in coordinates of certain types of beads from a lammps trajectory file (slightly modified version of original ```pppmd.read_lammpstrj``` function. <a name="special_read"/>

    **Parameters**: ```fname```, ```types```, ```num_frames```, ```skip_beginning```, ```skip_between```
    - ```fname```: filename string or 'stdin' (or a value that evaluates to false) for reading from standard in 
    - ```types```: *list* of atom types to extract coordinates for
    - ```num_frames```: optional number of frames to read before stopping, defaults to reading in all frames
    - ```skip_beginning```: skip this many frames at the beginning of the dump file
    - ```skip_between```: skip this many frames between saved frames

    **Returns**: ```r```, ```ir```, ```timestep```, ```box_bounds```, ```id2type```, ```id2index```
    - ```r```: num_frames by num_atoms+1 by 3 array of wrapped and unscaled coordinates (indexed by frame number then atom id)
    - ```ir```: num_frames by num_atoms+1 by 3 array of image flags
    - ```timestep```: num_frames length array of timesteps
    - ```box_bounds```: 3D array to store boundaries of the box, indexed by frame, x/y/z, then lower/upper
    - ```id2type```, ```id2index```: num_atoms+1 length arrays to map atom id to type and index id (if available, may be None)

    This function also outputs the amounts of total beads and selected beads as follows:
    ```
    ===============Summary===============
    Total number of beads =  num_atoms
    Total number of selected beads = num_types
    ```

2. ```wrap```: Simply wrap coordinates; this is needed because atoms can be slightly outside the periodic boundary even if they are dumped as wrapped coordinates from LAMMPS (see [source](https://lammps.sandia.gov/threads/msg32219.html)). <a name="wrap"/>

    **Parameters**: ```r```, ```box_bounds```
    - ```r```, ```box_bounds``` are outputs from ```special_read``` function

    **Returns**: ```r```, ```boxsize```, ```bound_lo```
    - ```r```: num_frames by num_atoms+1 by 3 of once-again wrapped coordinates
    - ```boxsize```: 3D array to store box dimensions *(does not account for changes in box size)*
    - ```bound_lo```: 3D array to store lower boundaries of the box *(does not account for changes in box size)*

3. ```binbox```: Bins the box and assigns each atom into the corresponding bin.  <a name="binbox"/>

    **Parameters**: ```r```, ```boxsize```, ```bound_lo```, ```dist_range```
    - ```r```, ```boxsize```, ```bound_lo``` are outputs from ```wrap``` function
    - ```dist_range```: bin size float number
    
    **Returns**: ```bin2id```, ```id2bin```, ```bins```
    - ```bin2id```: number of bins in x by number of bins in y by number of bins in z by num_atoms+1 by num_frames list to map bin to atom id
    - ```id2bin```: num_frames by num_atoms+1 by 3 array to map atom id to bin
    - ```bins```: 3D array to store number of bins in x/y/z

4. ```buildnlist```: Builds a neighbor list for each bead and return a id2neighbors: num_frames by num_atoms+1 by 3 array of neighbor list.   <a name="buildnlist"/>
    
    **Parameters**: ```r```, ```bin2id```, ```id2bin```, ```bins```, ```boxsize```, ```id2type```, ```dist_range```, ```nearest```
    - ```r```, ```bin2id```, ```id2bin```, ```bins```, ```boxsize```, ```id2type```, ```dist_range``` are outputs from functions above or already defined above
    - ```nearest```: Boolean to determine whether to have only the nearest bead in the neighbor list

    **Returns**: ```id2neighbors```
    - ```id2neighbors```: num_frames by num_atoms+1 by number of neighbors array to map atom id to neighbor list
    
5. ```ipcorr```: Calculate ion pair autocorrelation function. The definition can be found in our [paper](https://pubs.acs.org/doi/abs/10.1021/acs.macromol.0c00216)  <a name="ipcorr"/>
    
    **Parameters**: ```r```, ```id2neighbor```
    - ```r```, ```id2neighbor``` are outputs from functions above or already defined above

    **Returns**: ```id2neighbors```
    - ```corr```: ion pair autocorrelation function versus timestep

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

*Example usages of kshen/ion-dynamics functions are available in [kshen/ion-dynamics/example](https://github.com/hall-polymers/pppmd2/tree/development/kshen/example).*
