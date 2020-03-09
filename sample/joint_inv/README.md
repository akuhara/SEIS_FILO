# Transdimensional joint inversion of surface wave and receiver function

## How to run

`mpirun -np 20 ../../bin/joint_inv joint_inv`
* The number after the "-np" flag represents the number of processes. Any numbers >= 1 is acceptable, but for small numbers, the convergence will be slow (i.e., requiring more numbers of iterations).

## Input files

* [joint_inv.in](https://github.com/akuhara/SEIS_FILO/blob/master/sample/joint_inv/joint_inv.in): the main parameter file
* [disper.in](https://github.com/akuhara/SEIS_FILO/blob/master/sample/joint_inv/disper.in): observation summary file for dispersion curves
* [rayleigh.0th](https://github.com/akuhara/SEIS_FILO/blob/master/sample/joint_inv/rayleigh.0th): observed data (phase & group velocities of the fundamental Rayleigh mode)
* `recv_func.in`: observation summary file for receiver functions
* `recv_func.sac`: observed data (P receiver function)

## How to make plots

1. `python ../../util/plot_recv_func.py joint_inv.in 1`, which produces `disper01.png`.

2. `python ../../util/plot_disper.py joint_inv.in 1`, which produces `recv_func01.png`.
