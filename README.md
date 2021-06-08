# SEIS_FILO (Beta version) test test

__SEISmological transdimensional inversion tools for Flat and Isotropic Layered structures in the Ocean__ 

[![Build Status](https://travis-ci.org/akuhara/SEIS_FILO.svg?branch=master)](https://travis-ci.org/akuhara/SEIS_FILO)
[![codecov](https://codecov.io/gh/akuhara/SEIS_FILO/branch/master/graph/badge.svg)](https://codecov.io/gh/akuhara/SEIS_FILO)
![GitHub](https://img.shields.io/github/license/akuhara/SEIS_FILO)
![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/akuhara/seis-filo)
[![Documentation Status](https://readthedocs.org/projects/seis-filo/badge/?version=latest)](https://seis-filo.readthedocs.io/en/latest/?badge=latest)

Copyright (C) 2019-2021 Takeshi Akuhara

![LOGO](./img/SEIS_FILO_LOGO.png)
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
    * Flat or spherical Earth model
    * Solve for noise standard deviation
    * Prior constraint applied or not
* __Easy to visualize__
    * Plot utilities also available

---

## Requirements
* [FFTW library](http://fftw.org/)
* [LAPACK library](http://www.netlib.org/lapack/)
* [Open MPI](https://www.open-mpi.org/)
* Some Python modules listed in [requirements.txt](https://github.com/akuhara/SEIS_FILO/blob/master/requirements.txt)

  
## Install
Type `make` in the `src` directory. Please edit the [Makefile](https://github.com/akuhara/SEIS_FILO/tree/master/src/Makefile) in accordance with your environment (i.e., compiler type and libarary paths). 

Alternatively, you can choose to use [docker container](https://hub.docker.com/r/akuhara/seis-filo) in which ready-to-use executable files are stored. 

## Quick Start with Sample
* [__disper_fwd__](https://github.com/akuhara/SEIS_FILO/tree/master/sample/disper_fwd): Surface wave forward computation (phase & group velocities)
* [__recv_func_fwd__](https://github.com/akuhara/SEIS_FILO/tree/master/sample/recv_func_fwd): Receiver function forward computation
* [__joint_inv__](https://github.com/akuhara/SEIS_FILO/tree/master/sample/joint_inv): Surface wave and receiver function joint inversion by RJMCMC

---

## For More Details
See [online documentation](https://seis-filo.readthedocs.io).



