# SEIS_FILO
SEISmological tools for Flat Isotropic Layered structure in the Ocean

Copyright (C) 2019 Takeshi Akuhara

---

## Programs included in this package
* Rayleigh wave forward computation (rayleigh_fwd)
* Rayleigh wave transdiemsnional inversion (rayleigh_inv)

## Install
Type `make` in the `src` directory.

---

## Rayleigh wave forward computation
### USAGE
`rayleigh_fwd [parameter file]`
### Parameter file
* Need to specify the following parameters in this file.
  * fmin, fmax, df (Minimum, maximum, and interval of frequency at which dispersion curve is calculated)
  * cmin, cmax, dc (Minimum, maximum, and interval of phase velocity used for root search)
  * vmod_in (file name for input velocity model)
  * ray_out (file name for output dispersion curve)
* Comment out by "#" 
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

## Rayleigh wave inversion

in preparation

---


