# SEIS_FILO (Alpha version)
SEISmological inversion tools for Flat and Isotropic Layered structures in the Ocean

# NOTE:
Currently, the software undergoes alpha testing: part of functions seem to work but not thoroughly tested yet. 

Copyright (C) 2019 Takeshi Akuhara

![plot](./img/plot.png)
---

## Quick Start

* __disper_fwd__: Surface wave forward computation (phase & group velocities)
* __recv_func_fwd__: Receiver function forward computation
* [__joint_inv__](https://github.com/akuhara/SEIS_FILO/tree/master/sample/joint_inv): Surface wave and receiver function joint inversion by RJMCMC
  
## Install
Type `make` in the `src` directory.

## Requirements
* [FFTW library](http://fftw.org/)
* [LAPACK library](http://www.netlib.org/lapack/)
* [Open MPI](https://www.open-mpi.org/)
---

## How to Use
See [Wiki](https://github.com/akuhara/SEIS_FILO/wiki).
