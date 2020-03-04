# Sample directory for joint inversion of P receiver functions and Rayleigh waves (phase & group velocities of the fundamental mode)

## How to run

`mpirun -np 20 ../../bin/joint_inv joint_inv`
* The number after the "-np" flag represents the number of processes. Any numbers >= 1 is acceptable, but for small numbers, the convergence will be slow (i.e., requiring more numbers of iterations).


## How to make plots

1. `python ../../util/plot_recv_func.py joint_inv.in 1`, which produces `disper01.png`.


2. `python ../../util/plot_disper.py joint_inv.in 1`, which produces `recv_func01.png`.