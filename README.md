# SEIS_FILO
SEISmological tools for Flat Isotropic Layered structure in the Ocean

## Rayleigh wave forward computation
### USAGE
`rayleigh_fwd [parameter file]`
### Parameter file

```
# <- can be used for comment out
# No space is allowed at both sides of "="
fmin=0.01 fmax=1.0  df=0.1  # Minimum, maximum, and interval of frequency at which dispersion curve is calculated
cmin=1.0  cmax=2.0  dc=0.01 # Minimum, maximum, and interval of phase velocity used for root search 
vmod_in=vmod.in             # file name for input velocity model
ray_out=ray.out             # file name for output dispersion curve
```

### Velocity model file (vmod_in)

```
3                 (Number of layers)
2.55 1.50 1.0 2.0 (Vp, Vs, Density, and Thickness of 1st layer)
3.40 2.00 1.0 2.0 (                                  2nd layer)
5.10 3.00 1.0 2.0 (                                  3rd layer)

```
---

## Install
Type `make` in the `src` directory.
