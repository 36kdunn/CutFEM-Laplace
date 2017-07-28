# CutFEM-Laplace
Solve -Lap(u)=f by running run_intExtFEM.m which calls intExtFEM.m.
Points on the interface Gamma are given in G, which is defined in run_intExtFEM.m.
In the current configuration, Dirichlet BCs are prescribed on the boundary of the square
and a zero jump condition is imposed on the interface.
The jump of the normal derivative is also set to be zero at the interface.
These boundary conditions can be changed by editing the corresponding files in the BCs folder.

All files in each folder are required to run the code.
Running run_intExtFEM.m in Matlab will add the necessary included folders to the default path.

