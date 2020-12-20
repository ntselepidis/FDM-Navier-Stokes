# Navier-Stokes Simulation

The code is divided in the following files:
```
├── generate_plot.m  <-- MATLAB script used to generate plots
├── main.f90         <-- Main program implementing Navier-Stokes simulation
├── Makefile
├── m_bc.f90         <-- Boundary conditions (used in m_sim.f90 and m_mg.f90)
├── m_fdm.f90        <-- Finite difference method (diffusion and advection computation)
├── m_mg.f90         <-- Multigrid solver
├── m_sim.f90        <-- Simulation module
├── m_stream.f90     <-- Stream function and velocity computation
├── parameters.txt   <-- Main parameter file
└── README.md
```
