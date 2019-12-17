# SEIS_FILO
SEISmological tools for Flat Isotropic Layered structure in the Ocean

Copyright (C) 2019 Takeshi Akuhara

---

## Programs included in this package
* __rayleigh_fwd__: Rayleigh wave forward computation
* __rayleigh_inv__: Rayleigh wave transdimensional inversion by RJMCMC

## Install
Type `make` in the `src` directory.

---

## Rayleigh wave forward computation
### USAGE
`rayleigh_fwd [parameter file]`
### Parameter file
* Need to specify the following parameters in this file.
  * __fmin__, __fmax__, __df__ (Minimum, maximum, and interval of frequency at which dispersion curve is calculated)
  * __cmin__, __cmax__, __dc__ (Minimum, maximum, and interval of phase velocity used for root search)
  * __vmod_in__ (file name for input velocity model)
  * __ray_out__ (file name for output dispersion curve)
* Comment out by "#" works fine.
```
# Parameter file example
fmin = 0.01 # Can add comment
fmax = 1.0  
df = 0.1  
cmin = 1.0  
cmax = 2.0  
dc = 0.01 
vmod_in = vmod.in            
ray_out = ray.out        
```

### Velocity model file (vmod_in)

* The first line should be the number of layers
* Vp, Vs, density and thickness of each layer should be listed in the successive lines
* Set Vs < 0 for the ocean layer (for now, only the topmost layer is allowed for this)
* Any comment out does not work for this file 
```
3                 
2.55 1.50 1.0 2.0 
3.40 2.00 1.0 2.0 
5.10 3.00 1.0 2.0 
```

### Dispersion curve file (ray_out)

* Each line contains frequency, phase velocity, and group velocity.
```
    0.0100    2.7122    2.6735
    0.1100    2.3291    1.9678
    0.2100    1.8198    1.1672
    0.3100    1.5359    1.1804
    0.4100    1.4420    1.2438
    0.5100    1.4049    1.2950
    0.6100    1.3890    1.3288
    0.7100    1.3819    1.3491
    0.8100    1.3785    1.3607
    0.9100    1.3770    1.3674
    1.0100    1.3762    1.3711
```
---

## Rayleigh wave inversion

### USAGE
`mpirun -np 20 rayleigh_inv [parameter file]`
* The numeric following '-np' indicates the number of processes for parallel computing.

### Parameter file
* Need to specify the following parameters in this file.
  * __n_iter__, __n_corr__, __n_burn__ (# of total iterations, iterations for sampling interval, and iterations in burn-in phase)
  * __dev_vs__, __dev_vp__, __dev_vz__ (Standard deviation for random walk along Vs, Vp and Z axes)
  * __vs_min__, __vs_max__, __vp_min__, __vp_max__, __z_min__, __z_max__ (Prior bounds, where uniform distribution is assumed)

* Comment out by "#" works fine.


```
# Parameter file example

# MCMC iteration
n_iter=30000
n_corr=100
n_burn=10000

# MCMC proposal
dev_vs=0.1
dev_vp=0.1
dev_z=0.1

# MCMC prior
vs_min=2.5
vs_max=5.0
vp_min=5.0
vp_max=8.5
z_min=0.0
z_max=30.0

# Parallel temerping
n_chain=5
n_cool=1
temp_high=10.0

# Transdiensional model space
k_min=1
k_max=21 

# Random number seeds
i_seed1=14222  
i_seed2=144444 
i_seed3=98767889 
i_seed4=22405559

# Observation file
obs_in=obs.in

# Rayleigh dispersion computation
cmin=2.5
cmax=4.5
dc=0.01

# Inversion setting
solve_vp=.true.
ocean_flag=.false.
ocean_thick=1.d0

# Output binning
nbin_z =50
nbin_vp = 25 
nbin_vs = 25 

```
---


