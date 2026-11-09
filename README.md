# gPLUTO
Code version: 0.91  (June 2026)

gPLUTO is a freely-distributed software for the numerical solution of mixed hyperbolic/parabolic systems of partial differential equations (conservation laws) targeting high Mach number flows in astrophysical plasma dynamics. 
The code is designed with a modular and flexible structure whereby different numerical algorithms can be separately combined to solve systems of conservation laws using the finite volume or finite difference approach based on Godunov-type schemes.
The present project replaces the legacy software PLUTO by providing Nvidia GPU support.

The code is written in the C++ programming language.

gPLUTO is a highly portable software and can run from a single workstation up to several thousands processors using the Message Passing Interface (MPI/NCCL) to achieve highly scalable parallel performance.

The software is developed at the [Dipartimento di Fisica, Torino University](https://www.unito.it/) in a joint collaboration with [INAF](https://www.oato.inaf.it/?lang=en), Osservatorio Astronomico di Torino and currently under the support of the  [SPACE Center of Excellence](https://www.space-coe.eu/).

## System requirements

Depending on the desired configuration (serial/parallel) different system requisites must be satisfied in order to compile and run gPLUTO. These consist of a minimal set plus additional libraries:

- Minimal prerequisites: C++ Compiler, GNU make, Python3
- Parallel builds: MPI library 


If you have a laptop with an NVIDIA GPU, you may want to use the GPU support by downloading and installing the NVIDIA compiler HPC toolkit:

[https://developer.nvidia.com/hpc-sdk]()

## Installation

Clone the git repository via SSH
```bash
git clone git@gitlab.com:PLUTO-code/gPLUTO.git
```
or HTTPS
```bash
git clone https://gitlab.com/PLUTO-code/gPLUTO.git
```
add the environment variable gPLUTO_DIR to your ./bashrc (or similar)

```bash
echo "export gPLUTO_DIR=/path/to/dir/gPLUTO" >> ~/.bashrc
source ~/.bashrc
```

## Quick start
gPLUTO can be quickly configured to run one of the several test problems provided with the distribution. Assuming that your system satisfies all the requirements described in the Requirements chapter (i.e. C++ compiler, Python, etc..) you can quickly setup gPLUTO in the following way:
1. Change directory to any of the test problems
```bash
cd $gPLUTO_DIR/Test_Problems/HD/Sod
```

2. Copy the header and initialization files from a configuration of our choice (e.g. #01):
```bash
cp definitions_01.hpp definitions.hpp
cp pluto_01.ini pluto.ini
```
3. Run the Python script using
```bash
python $gPLUTO_DIR/setup.py
```
and select "Setup problem" from the main menu. You can choose (by pressing Enter) or modify the default setting using the arrow keys.

4. Once you return to the main menu, select "Change makefile", choose a suitable makefile (e.g. linux.g++.defs) and press enter. All the information relevant to the specific problem should now be stored in the four files init.cpp (assigns initial condition and user-supplied boundary conditions), pluto.ini (sets the number of grid zones, Riemann solver, output frequency, etc.), definitions.hpp (specifies the geometry, number of dimensions, reconstruction, time stepping scheme, and so forth) and the makefile.

5. Exit from the main menu ("Quit" or press 'q') and type
```bash
make -j
```
to compile the code.

6. You can now run the code by typing
```bash
./pluto
```
At this point, gPLUTO reads the initialization file pluto.ini and starts integrating. The run should take a few seconds (or less) and the integration log should be dumped to screen.

## Code Overview
### Physics modules
At present, gPLUTO can solve different systems of conservation laws in 1, 2 or 3 dimensions using Cartesian coordinates. The four basic physics modules are

- Classical hydrodynamics (Euler or Navier-Stokes equations, HD);
- Magnetohydrodynamics (MHD);
- Special Relativistic hydrodynamics (RHD);
- Special Relativistic MHD (RMHD);

On top of these, different equations of state are available
- Isothermal
- Ideal
- Synge

as well as fixed gravity.


Equations are discretized and solved on a structured static mesh.


### Available numerical methods.

Modularity applies not only to the type of equations being solved but also to the capability of enabling different numerical algorithms to be employed and combined in different contexts and different components. 

gPLUTO solves the fluid equations by means of shock-capturing Godunov-type methods adopting a conservative discretization based on finite volume or finite difference methods. 
In this formulation interface fluxes are computed by solving a Riemann problem between left and right interface states. 
A typical time step cycle consists of a
1) an explicit time-marching algorithm
2) a piece-wise reconstruction scheme
3) a Riemann solver.


MHD and relativistic MHD equations may be solved by selecting different strategies to enforce the divergence-free condition, including Constrained transport (CT) / hyperbolic divergence cleaning / Powell's 8 wave formulation.

Please refer to the documentation page for more information.


## Terms of use.

gPLUTO is distributed freely under the GNU general public license. Codes development and support requires a great deal of work and for this reason we expect gPLUTO to be referenced and acknowledged by authors who use it for their publications. Co-authorship may be solicited for those publications demanding considerable additional support and/or changes to the code.



## Method Papers

[4th-Order Finite Volume](https://ui.adsabs.harvard.edu/abs/2024JCoPh.49912701B/abstract) method paper: "A 4th-order accurate finite volume method for ideal classical and special relativistic MHD based on pointwise reconstructions", Berta et al., JCoPh (2024) 49912701B

[Upwind Constrained Transport](https://ui.adsabs.harvard.edu/abs/2021JCoPh.42409748M/abstract) method paper: "Systematic construction of upwind constrained transport schemes for MHD",  A. Mignone, L. Del Zanna, JCoPh (2021) 42409748M

[Resistive Relativistic MHD](https://ui.adsabs.harvard.edu/abs/2019MNRAS.486.4252M/abstract) method paper: "A constrained transport method for the solution of the resistive relativistic MHD equations",  Mignone et al., MNRAS (2019) 486, 4252M

[Lagrangian Particles](https://ui.adsabs.harvard.edu/abs/2018ApJ...865..144V/abstract) method paper: "A Particle Module for the PLUTO Code. II. Hybrid Framework for Modeling Nonthermal Emission from Relativistic Magnetized Flows", Vaidya et al., ApJ (2018) 865, 144V


[Finite Difference](https://ui.adsabs.harvard.edu/abs/2010JCoPh.229.5896M/abstract) method paper: "High-order conservative finite difference GLM-MHD schemes for cell-centered MHD", Mignone et al., JCoPh (2010) 229, 5896M


[Static grid](https://ui.adsabs.harvard.edu/abs/2007ApJS..170..228M/abstract) (original) method paper: "PLUTO: A Numerical Code for Computational Astrophysics", Mignone et al., ApJS (2007) 170, Issue 1, pp. 228-242
