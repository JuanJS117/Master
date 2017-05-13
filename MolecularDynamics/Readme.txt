MOLECULAR DYNAMICS OF A LENNARD JONES FLUID -- By Juan Jimenez Sanchez


PURPOSE

This is the final work of the 'Computational Methods in Condensed Matter Physics and
Biomolecules' course done at 'Master in Condensed Matter and Biological Systems Physics'
in UAM, Madrid.

The aim is to model the behaviour of a Lennard-Jones fluid using Molecular Dynamics,
in order to study some hydrodynamical properties at nanoscale, by checking Time
Autocorrelation Functions. Three transport coefficients are studied here:

    1. Diffusion coefficient -> It comes from integrating the Velocity Autocorrelation
       Function (VACF) over time, and multiplying it by a 1/3 factor.

    2. Friction coefficient -> It comes from integrating the Force Autocorrelation Function
       (FACF) over time, and multiplying it by a 1/3T factor.

    3. Kinematic viscosity -> It is not so simple as the previous coefficients. We first need
       to get the momentum field autocorrelation function (GACF) in Fourier space. Then,
       following this relation:

       <g(k,t).g(k,0)> = <g(k,0).g(k,0)>exp(-visc.k².t) --> where visc is the kinematic VISCOSITY

       We need to fit the ln of the normalized GACF to get the kinematic viscosity. k is related
       to the wavenumber of our system, so the kinematic viscosity will be related to k:

       visc_k = visc_0 - A.k²

       The true kinematic viscosity is visc_0. To get it, we need to fit the visc values obtained
       to the different k used to calculate GACF.

In this work, we checked that the Einstein-Smoluchowski relation is accomplished for the calculations
done. This sets this approach as a good one, as long as all results are consistent with evidence.


FILES

In this repository we will find:

1. 'md_NVE.c' -> The main code. A program written in C which runs a Molecular Dynamics
    simulation in NVE. To keep temperature constant around a fixed value, the initial coordinates
    are read from a file containing an equilibrium Montecarlo configuration at the desired temperature.
    If keeping temperature constant is not a must, an equally-spaced lattice of particles may also be
    created. Also, the initial velocities are randomly distributed following a Boltzmann distribution
    at the desired temperature. To obtain such uniformly distributed random values, the polar Box-Muller
    algorithm is implemented. The position, velocity and force integration is done using the Velocity
    Verlet algorithm. Each N (N = NSTEP/NSAMP) simulation steps, the program calculates
    the potential, kinetic and total energy, the system temperature and pressure, and the
    velocity, force, and momentum field autocorrelation functions. It also prints the current
    system state on the screen. In the end, the velocity and force autocorrelation functions
    are integrated over time, in order to get the diffusion and the friction coefficients
    (See Green-Kubo relations and Transport Coefficients as integrals of Time Autocorrelation
    Functions). Additionally, the program stores each N steps all the positions, velocities and
    forces for each particle in a file. All those trajectory files may be used later to make
    a gif of the full Molecular Dynamics run.

    To run the program, please take a look first at the simulation parameters, and change them as
    you wish, but keeping in mind that the actual values ensure stability. In Linux, the shell
    command to compile it is:

    gcc md_NVE.c -lm

    This will create the output file 'a.out'. To run it and store system info, run the command:

    ./a.out > file.out &

    For a 200 particle system, the simulation will take around 8-11 minutes.

2. 'plot_trajectories.R' -> The main analysis code. A program written in R to check every
    variable of the Molecular Dynamics simulation. In the end, the trajectory files are all
    read and plotted in a 3D scatterplot. An image is stored. All those images may be joined
    together into a GIF with the following shell command:

    convert *.png -delay 3 -loop 0 md.gif

    This program is full-purpose, but more detailed codes have been also written in order to
    make clear and understandable plots.
    Rstudio is recommended to quickly access every plot. However, you can also call the
    program via shell using the following command:

    ./plot_trajectories.R

    Simply take care of the simulation output files, in order to keep the same filenames in
    this program.

3. 'MD_test_plots.R' -> This program written in R is intended to check how well a Molecular
    Dynamics simulation conserves the temperature, the energy, and the momentum along each
    direction. It can read from 1 to 5 different simulation files.

4. 'Einstein_Stokes_Relations.R' -> This program, also written in R, is able to read the
    FACF and VACF files of 1-5 different files, and check whether the Einstein-Smoluchowski
    and the Stokes relations are accomplished or not.

5. 'get_kinvisc.R' -> This program in R is intended to process the info contained in the GACF
    file, and retrieve the kinematic viscosity. The way to get it is described before; in this
    code all I do is to implement the explained scheme.

6. 'md.gif' -> An already done GIF of the first successful Molecular Dynamics simulation.
    It shows 200 particles moving along the box, with an initial configuration consisting
    in an equally-spaced 3D lattice of particles.

7. '*.png' -> All images are plots obtained from the R codes. The filename reveals its content.


IMPROVEMENTS

Of course, there is much more to do in order to improve the codes. None has been designed worrying
about user interface, as it was not intended for the course. But all of them are well commented,
so anyone interested in using them will find easy to change system parameters.

Additionally, the runtime may be largely improved by implementing a neighbour search subroutine.
But the runtime for 200 particles systems are not very high (8-11 min), and it is my first
code done in C; both the complexity and the available time took me to avoiding the neighbour
search subroutine. In further improvements a linked cell list subroutine may be written.

For any question, please contact me at my institutional e-mail:
juan.jimenezs@estudiante.uam.es
