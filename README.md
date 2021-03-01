# SCFT_CUDA
This is an SCFT program implemented by the CUDA technique.

# Purpose
This program is to calculate the free energy and density profile of a diblock copolymer melt in equilibrium (canonical ensemble). The mathematical description of the program is written in “Stable and Unstable Phases of a Diblock Copolymer Melt.”
Matsen, M. W., and M. Schick. Physical Review Letters, vol. 72, no. 16, 1994, pp. 2660–2663.

# Usage
1. Set the configuration.
The configuration of SCFT locates at ab.dat, and is read by the source file /scr/init.cpp;

intag=1024    // intag is a marker that tells the program how to initialize the conjugation field wA and wB. Here 1024 means the wA and wB is initialized from file /Density/phi_1.dat
xN=30.0 // \chi N is the incompatible parameter between A block and B block.
fa=0.78  // fa is the volume fraction of A block and the B block is set automatically as 1-fa.
lx=29.30, ly=14.5,  lz=2.4494897427831781 //  The three dimensions of the box in the unit of Rg.
Nx=128, Ny=64, Nz=1  // The discretization of the box in three dimensions
ds0=0.01   // The discretization of the polymer chain, the number of segment per chain is 1/ds0. For example, here the number of segments per chain is 100.
wopt=0.2   // Normally, you can ignore this parameter. It is the mixing parameter for the self-consistent iteration. When the program failed to converge, you may try to decrease this parameter.
ErrorinCompMax=0.0000001 // The Maximal value of the incompressibility that the program can tolerate. Usually, it is set to 10^-7.
batch=1 // This program can deal with multiple calculations at the same time, the batch gives the number.

2. Type make to compile to program.

3. Type ./scft -Steps_N=10000. The SCFT iteratively runs 10000 steps or stops when the incompressibility is smaller than ErrorinCompMax=0.0000001. 
you may also type ./scft, the program runs 1 step by default.

4. The screen outputs the following information:
The iterative steps; free energy; the interfacial energy; the conformational entropy; the incompressibility error;
The density profile locates at /Density/pha_1.dat. You can visulize it by paraview, i.e.,
a. Type gcc vtk.c in the command line.
b. Type ./a.out pha_1.dat.
c. Copy or drag the generated file out.vtk to the paraview window.
d. Press the apply button.

The program outputs every 100 steps by default, you can change the frequency by adding flag: -Check_T=#steps.


