## PPPMD dump_tool package
**Developed by Nicholas Liesen**

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
