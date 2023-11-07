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
n_iter           = 200  # Number of total iterations
n_burn           = 100  # Number of iterations for burn-in period
n_corr           = 1	   # Interval of iterations to evaluate models

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
rayleigh.0th
R 0 period
30 1.0 1.0
0.3 4.5 0.05
0.005 0.09 0.03
0.005 0.1 0.01
0.005 0.09 0.03
0.005 0.09 0.03
#------------------------------------------------
# If multiple inputs (n_disp >= 2), write them below.
#------------------------------------------------
#rayleigh.1th
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
recv_func.sac
8.0 0.05 0.05
P
0.0 2.5
0.005 0.05 0.005
T T 0.001
#----------------------------------------------------------
# If multiple inputs (n_rf <= 2), write them below.
#----------------------------------------------------------
#recv_func2.sac
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

1.5013 T 1.4521 F 0.4141 F # 1.0000
1.5506 T 1.4220 F 0.4077 F # 2.0000
1.7351 T 1.2729 F 0.4306 F # 3.0000
1.9642 T 1.1640 F 0.4671 F # 4.0000
2.3346 T 1.1823 F 0.4452 F # 5.0000
2.8300 T 1.5285 F 0.5179 F # 6.0000
3.1888 T 2.0098 F 0.5738 F # 7.0000
3.3816 T 2.4087 F 0.6466 F # 8.0000
3.5165 T 2.6757 F 0.6860 F # 9.0000
3.6116 T 2.8739 F 0.7789 F # 10.0000
3.6700 T 3.0576 F 0.7876 F # 11.0000
3.7606 T 3.1680 F 0.8124 F # 12.0000
3.7919 T 3.3388 F 0.8214 F # 13.0000
3.8527 T 3.4106 F 0.8734 F # 14.0000
3.9082 T 3.4702 F 0.8716 F # 15.0000
3.9215 T 3.5385 F 0.8878 F # 16.0000
3.9424 T 3.5949 F 0.8676 F # 17.0000
3.9773 T 3.6727 F 0.9240 F # 18.0000
3.9495 T 3.6801 F 0.8705 F # 19.0000
3.9735 T 3.7078 F 0.9387 F # 20.0000
3.9802 T 3.7549 F 0.9093 F # 21.0000
4.0058 T 3.8055 F 0.9593 F # 22.0000
3.9953 T 3.8220 F 0.8982 F # 23.0000
4.0158 T 3.8250 F 0.9131 F # 24.0000
4.0536 T 3.8668 F 0.9106 F # 25.0000
4.0429 T 3.8985 F 0.9289 F # 26.0000
4.0266 T 3.8743 F 0.9052 F # 27.0000
4.0321 T 3.8827 F 0.9115 F # 28.0000
4.0409 T 3.8956 F 0.8837 F # 29.0000
4.0822 T 3.9208 F 0.8999 F # 30.0000

```

### Receiver function data file
* The filename must be specified in the receiver function data summary file.
* Use SAC file format.

## How to make plots

1. `python ../../../util/plot_disper.py joint_inv.in 1`, which produces `disper01.png`.

2. `python ../../../util/plot_recv_func.py joint_inv.in 1`, which produces `recv_func01.png`.
