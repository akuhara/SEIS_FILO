# SEIS_FILO (Alpha version)
__SEISmological inversion tools for Flat and Isotropic Layered structures in the Ocean__ 

Copyright (C) 2019 Takeshi Akuhara

---

The SEIS_FILO program package aims to carry out transdimensional joint inversion of surface waves and receiver functions for ocean-bottom observatories. The main features are: 

* __Transdimensional MCMC__
    * Automatically determines the number of unknowns (i.e., the number of layers in the structure model).
    * Uncertainty estimates
* __Parallel computing__
    * Paralle tempring technique is available  
* __Ocean layer__
    * Acoustic solution for the ocean layer and elastic solution for solid layers
* __Multiple inputs__
    * Dispersion curves of the fundamental and higher mode Rayleigh waves (Love wave is not implimented yet but is planned)
    * Rayleigh wave ellipticity & admittance (planned)
    * P and S receiver functions with different filter parameters and ray parameters
    * Works fine with any combination of inputs above 
* __Flexible inversion setting__
    * Vp solved or fixed
    * Ocean exists or not
    * Flat or spherical model (planned)
    * Solve for noise standard deviation (planned)
* __Easy to visualize__
    * Plot utilities also available
    
__NOTE: Currently, the software undergoes alpha testing: part of functions seem to work but not thoroughly tested yet. Anythings in codes and documents can change without notification.__


---

## Requirements
* [FFTW library](http://fftw.org/)
* [LAPACK library](http://www.netlib.org/lapack/)
* [Open MPI](https://www.open-mpi.org/)

  
## Install
Type `make` in the `src` directory. Please edit the `Makefile` in accordance with your environment (i.e., compiler type and libarary paths). 


## Quick Start
* __disper_fwd__: Surface wave forward computation (phase & group velocities)
* __recv_func_fwd__: Receiver function forward computation
* [__joint_inv__](https://github.com/akuhara/SEIS_FILO/tree/master/sample/joint_inv): Surface wave and receiver function joint inversion by RJMCMC

---

## How to Use
See [Wiki](https://github.com/akuhara/SEIS_FILO/wiki).

![plot](./img/plot.png)
