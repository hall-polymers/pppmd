## Original PPPMD Package developed by Jonathan R. Brown

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
