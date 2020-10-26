# Requirements

The following libraries must be installed. 

* [FFTW library](http://fftw.org/)
* [LAPACK library](http://www.netlib.org/lapack/)

The appropriate locations of the above libraries must be specified in `src/Makefile`. 

* [Open MPI](https://www.open-mpi.org/)

The MPI compiler must be linked to the GNU Fortran compiler (_i.e._, `gfortran`). The programs have been confirmed to work using the GNU Fortran compiler version 4.8.5. The older version may fail to compile the programs because the codes are witten with modern Fortran features.

Plot utilities require Python3 and popular modules including numpy, matplotlib, pandas, and seaborn. See `requirements.txt` for a full list of modules required.
