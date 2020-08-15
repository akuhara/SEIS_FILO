# Transdimensional inversion of receiver function with prior constraint

## How to run

`mpirun -np 20 ../../../bin/joint_inv joint_inv.in`
* The number after the "-np" flag represents the number of processes. Any numbers >= 1 is acceptable, but for small numbers, the convergence will be slow (i.e., requiring more numbers of iterations).

## Input files

* [joint_inv.in](https://github.com/akuhara/SEIS_FILO/blob/master/sample/joint_inv/case_2_recv_func_with_ref_vmod/joint_inv.in): main parameter file
* [recv_func.in](https://github.com/akuhara/SEIS_FILO/blob/master/sample/joint_inv/case_2_recv_func_with_ref_vmod/recv_func.in): data summary file for receiver functions
* [recv_func.sac](https://github.com/akuhara/SEIS_FILO/blob/master/sample/joint_inv/case_2_recv_func_with_ref_vmod/recv_func.sac): data file (P receiver function)

## How to make plots

1. `python ../../../util/plot_recv_func.py joint_inv.in 1`, which produces `disper01.png`.

