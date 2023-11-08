# Transdimensional Joint Inversion 

## How to run

`mpirun -np 20 ../../bin/joint_inv joint_inv.in`
* The number after the "-np" flag represents the number of processes.

## Inputs

### Main parameter file

* The filename must be specified as the first command line argument.
* The following can be used as a template.

```
#-----------------------------------------------------------------------
# Sample int to joint_inv of SEIS_FILO
#-----------------------------------------------------------------------

#--- Required parameters -----------------------------------------------

# Input data summary files
recv_func_in     = recv_func.in 
disper_in        = disper.in

# Model setting
solve_anomaly    = F       # Solve for velocity anomaly (T) or absolute velocity (F)
ref_vmod_in      =         # Filename of reference velocity model
solve_vp         = F       # Vp is solved (T) or not (F)
solve_rf_sig     = T       # Receiver function error is solved (T) or not (F)
solve_disper_sig = T       # Dispersion curve error is solved (T) or not (F)
is_sphere        = F       # Earth flattening is applied (T) or not (F)
is_ocean         = T       # Ocean layer exists (T) or not (F)
ocean_thick      = 2.0     # Ocean layer thickness [km]
vp_bottom        = 8.1     # Bottom half-space Vp [km/s]
vs_bottom        = 4.7     # Bottom half-space Vs [km/s]
rho_bottom       = 3.4     # Bottom half-space density [g/cm^3]

# Iteration numbers
n_iter           = 2000000  # Number of total iterations
n_burn           = 1000000  # Number of iterations for burn-in period
n_corr           = 1000	   # Interval of iterations to evaluate models

# Parallel temerping
n_chain        = 5         # Number of MCMC chains per process 
n_cool         = 1         # Number of non-tempered MCMC chains per process
temp_high      = 30.0      # The hottest temperature

# Prior probability 
k_min = 1                  # Lower limit of layer number
k_max = 21                 # Upper limit of layer number
# NOTE: k_min <= k < k_max
vs_min  = 1.0              # Lower limit of Vs [km/s]
vs_max  = 5.0              # Upper limit of Vs [km/s]
dvs_sig = 0.5              # Prior width for Vs anomaly [km/s]
vp_min  = 2.0              # Lower limit of Vp [km/s]
vp_max  = 8.5              # Upper limit of Vp [km/s]
dvp_sig = 0.5              # Prior width for Vp anomaly [km/s]
z_min   = 2.0              # Lower limit of layer interface depth [km]
z_max   = 23.0             # Upper limit of layer interface depth [km]

 
# MCMC proposal
dev_vs = 0.03              # Standard deviation for Vs random walk [km/s]
dev_vp = 0.08              # Standard deviation for Vp random walk [km/s]
dev_z  = 0.2               # Standard deviation for depth random walk [km]

# Output format
n_bin_vs  = 100            # Number of bin for Vs
n_bin_vp  = 100            # Number of bin for Vp
n_bin_z   = 100            # Number of bin for depth 
n_bin_sig = 100            # Number of bin for data error

#----- Optional parameters ----------------------------------------------

# For different planets
#r_earth = 6371.0           # Earth radius [km]

# Random number seeds
i_seed1 = 99999999  
i_seed2 = 33333333
i_seed3 = 24242424 
i_seed4 = 99887622

#-----------------------------------------------------------------------

```

### Dispersion curve data summary file

* The filename must be specified in the parameter file via a parameter "disper_in".
* The following can be used as a template.
* Parameters must apper in this order.

```
#------------------------------------------------
# Sample input to joint_inv of SEIS_FILO
#------------------------------------------------
# n_disp
1
#------------------------------------------------
# file name
# disper_phase, n_mode, freq_or_period
# nx, xmin, dx
# cmin, cmax, dc
# sig_c_min, sig_c_max, dev_sig_c
# sig_u_min, sig_u_max, dev_sig_u
# sig_hv_min, sig_hv_max, dev_sig_hv
# sig_ra_min, sig_ra_max, dev_sig_ra
#------------------------------------------------
disper_obs.txt
R 0 period
11 10.0 2.0
0.3 4.5 0.05
0.005 0.09 0.03
0.005 0.09 0.03
0.005 0.09 0.03
0.005 0.09 0.03
#------------------------------------------------
# If multiple inputs (n_disp >= 2), write them below.
#------------------------------------------------
#disper_obs2.txt
#R 1 period
#30 1.0 1.0
#0.3 4.5 0.05
#0.005 0.09 0.03
#0.005 0.1 0.01
#0.005 0.09 0.03
#0.005 0.09 0.03
```


### Receiver function data summary file
* The filename must be specified in the parameter file via a parameter "recv_func_in".
* The following can be used as a template.
* Parameters must apper in this order.

```
#----------------------------------------------------------
# Sample input to joint_inv of SEIS_FILO
#----------------------------------------------------------
# n_rf
1
#
#----------------------------------------------------------
# File name 
# a_gauss, rayp, delta
# rf_phase
# t_start, t_end
# sig_rf_min, sig_rf_max, dev_sig_rf
# deconv_flag, correct_amp, damp
#----------------------------------------------------------
recv_func_obs.sac
8.0 0.05 0.05
P
0.0 2.5
0.005 0.05 0.005
T T 0.001
#----------------------------------------------------------
# If multiple inputs (n_rf <= 2), write them below.
#----------------------------------------------------------
#recv_func_obs2.sac
#8.0 0.06 0.05
#P
#0.0 3.5
#0.005 0.05 0.005
#T T 0.001

```

### Dispersion curve data file
* The filename must be specified in the dipersion curve data summary file.
* The following can be used as a template

```
#-------------------------------------------------------------
# Sample Observation Data File:
# Input to joint_inv of SEIS_FILO package
#-------------------------------------------------------------

# Format
# 1st column: phase velocity (km/s)
# 2nd       : data availability (T or F)
# 3rd       : group velocity (km/s)
# 4th       : data availability (T or F)
# 5th       : H/V ratio
# 6th       : data availability (T or F)
# 7th       : Admittance (km/GPa)
# 8th       : data availability (T or F)

3.6100 T 2.8533 T 0.7668 T 0.9933 T # 10.0000
3.7823 T 3.1762 T 0.8556 T 1.5467 T # 12.0000
3.8545 T 3.4303 T 0.8396 T 2.2434 T # 14.0000
3.9076 T 3.5896 T 0.8578 T 2.9880 T # 16.0000
4.0121 T 3.6652 T 0.9221 T 3.8243 T # 18.0000
3.9564 T 3.7565 T 0.8941 T 4.7970 T # 20.0000
4.0117 T 3.7971 T 0.8721 T 5.8646 T # 22.0000
4.0221 T 3.8531 T 0.8810 T 7.0496 T # 24.0000
4.0270 T 3.8690 T 0.8896 T 8.2722 T # 26.0000
4.0280 T 3.9286 T 0.8809 T 9.6848 T # 28.0000
4.0618 T 3.9304 T 0.9222 T 11.1202 T # 30.0000


```

### Receiver function data file
* The filename must be specified in the receiver function data summary file.
* Use SAC file format.

## How to make plots

1. `python ../../../util/plot_disper.py joint_inv.in 1`, which produces `disper01.png`.

2. `python ../../../util/plot_recv_func.py joint_inv.in 1`, which produces `recv_func01.png`.
