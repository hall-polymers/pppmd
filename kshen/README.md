## PPPMD Ion Dynamics Package
**Developed by Kuan-Hsuan Kevin Shen**

This sub-package includes python analysis scripts for ion-containing systems, including ion pairing/clustering and ion conductivity analyses.

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