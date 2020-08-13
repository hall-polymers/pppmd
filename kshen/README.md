## PPPMD Ion Dynamics Package
**Developed by Kuan-Hsuan Kevin Shen**

This sub-package includes python analysis scripts for ion-containing systems, including ion pairing/clustering and ion conductivity analyses.

1. ```special_read```: Allows you to read in coordinates of certain types of beads from a lammps trajectory file (slightly modified version of original ```pppmd.read_lammpstrj``` function. This function also outputs the amounts of total beads and selected beads as follows
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
