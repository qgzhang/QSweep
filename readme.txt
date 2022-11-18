
--------------------------------------------------------------------------------------------------------
 This Demo is for the ICCV 2017 submission "Quasiconvex Plane Sweep for Triangulation with Outliers"
--------------------------------------------------------------------------------------------------------
1. The entry is 'main.m'.
2. Developed on Matlab 2016a, Windows 10 64-bits PC.


--------------------------------------------------------------------------------------------------------
Workflow of 'main.m' 
--------------------------------------------------------------------------------------------------------
1. Load the dataset 'Alcatraz Courtyard'.
2. Do triangulation using Polyhedron collapse and Q-sweep method.
3. Summarise run time and average converged error for each method.
4. Show 3D structrues generated respectively by Polyhedron collapse and Q-sweep.


--------------------------------------------------------------------------------------------------------
c-mex implementation of the function 'STEPSIZE' 
--------------------------------------------------------------------------------------------------------
1. 'c_stepsize.c' is for the 'STEPSIZE' function implemented in C. 
2. Three pre-compiled executables are included for Window 64-bit, Linux 64-bit and Mac 64-bit.