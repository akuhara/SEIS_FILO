# SEIS_FILO 

__SEISmological transdimensional inversion tools for Flat and Isotropic Layered structures in the Ocean__ 

[![Build Status](https://travis-ci.org/akuhara/SEIS_FILO.svg?branch=master)](https://travis-ci.org/akuhara/SEIS_FILO)
[![codecov](https://codecov.io/gh/akuhara/SEIS_FILO/branch/master/graph/badge.svg)](https://codecov.io/gh/akuhara/SEIS_FILO)
![GitHub](https://img.shields.io/github/license/akuhara/SEIS_FILO)
![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/akuhara/seis-filo)
[![Documentation Status](https://readthedocs.org/projects/seis-filo/badge/?version=latest)](https://seis-filo.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4082670.svg)](https://doi.org/10.5281/zenodo.4082670)

Copyright (C) 2019-2021 __Takeshi Akuhara__[![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0002-6129-8459)

___



The SEIS_FILO program package aims to carry out transdimensional joint inversion of surface waves and receiver functions for ocean-bottom observatories. The main features are: 

* __Transdimensional MCMC__
* __Parallel tempering__
* __Ocean layer__
* __Supported input types__
    * Dispersion curves of the fundamental and higher mode Rayleigh waves
    * Rayleigh wave ellipticity
    * P receiver functions
    * S receiver functions
* __Model parameters__
    * Absolute Vp & Vs
    * Vp and Vs anomalies relative to the reference
    * Number of layers
    * Layer depths
    * Standard deviation of noise

---

## Requirements
* [FFTW library](http://fftw.org/)
* [LAPACK library](http://www.netlib.org/lapack/)
* [Open MPI](https://www.open-mpi.org/)
* Some Python modules listed in [requirements.txt](https://github.com/akuhara/SEIS_FILO/blob/master/requirements.txt)

  
## Install
Type `make` in the `src` directory. Please edit the [Makefile](https://github.com/akuhara/SEIS_FILO/tree/master/src/Makefile) in accordance with your environment (i.e., compiler type and libarary paths). 

Alternatively, you can choose to use [docker container](https://hub.docker.com/r/akuhara/seis-filo) in which ready-to-use executable files are stored. 

## Quick Guidance
* [__disper_fwd__](https://github.com/akuhara/SEIS_FILO/tree/master/sample/disper_fwd): Surface wave forward computation (phase & group velocities)
* [__recv_func_fwd__](https://github.com/akuhara/SEIS_FILO/tree/master/sample/recv_func_fwd): Receiver function forward computation
* [__joint_inv__](https://github.com/akuhara/SEIS_FILO/tree/master/sample/joint_inv): Surface wave and receiver function joint inversion by RJMCMC
* Online documentation is also available [here](https://seis-filo.readthedocs.io).

## Publications

### Rayleigh wave dispersion curves
* Yamaya, L., Mochizuki, K., Akuhara, T., Nishida, K. (2021). Sedimentary structure derived from multi-mode ambient noise tomography with dense OBS network at the Japan Trench. _Journal of Geophysical Research: Solid Earth_, 126(6), e2021JB021789. https://doi.org/10.1029/2021JB021789

### S receiver functions
* Akuhara, T., Nakahigashi, K., Shinohara, M., Yamada, T., Shiobara, H., Yamashita, Y., et al. (2021). Lithosphere–asthenosphere boundary beneath the Sea of Japan from transdimensional inversion of S-receiver functions. Earth, Planets and Space, 73(1), 171. https://doi.org/10.1186/s40623-021-01501-5

___

![LOGO](./img/SEIS_FILO_LOGO.png)