# SEIS_FILO
SEISmological tools for Flat Isotropic Layered structure in the Ocean
Copyright (C) 2019 Takeshi Akuhara

## Rayleigh wave forward computation
### USAGE
`rayleigh_fwd [parameter file]`
### Parameter file
* Need to specify the following parameters in this file.
  * fmin, fmax, df
  * cmin, cmax, dc
  * vmod_in
  * ray_out

```
# <- can be used for comment out
# No space is allowed at both sides of "="
fmin=0.01 fmax=1.0  df=0.1  # Minimum, maximum, and interval of frequency at which dispersion curve is calculated
cmin=1.0  cmax=2.0  dc=0.01 # Minimum, maximum, and interval of phase velocity used for root search 
vmod_in=vmod.in             # file name for input velocity model
ray_out=ray.out             # file name for output dispersion curve
```

### Velocity model file (vmod_in)

* The first line should be the number of layers
* Vp, Vs, density and thickness of each layer should be listed in the successive lines
* Set Vs < 0 for the oceanic layer (for now, only the topmost layer is allowed for this)
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
   1.0000000000000000E-002   2.7121874999999998        2.6735455229920699     
  0.11000000000000000        2.3291406249999995        1.9678221055625724     
  0.21000000000000002        1.8198437499999993        1.1671994294820165     
  0.31000000000000005        1.5359374999999993        1.1804495815365621     
  0.41000000000000003        1.4420312499999990        1.2437814625397263     
  0.51000000000000001        1.4049218750000001        1.2950275266526103     
  0.61000000000000010        1.3889843750000004        1.3287697749037055     
  0.71000000000000008        1.3818749999999997        1.3490536150158645     
  0.81000000000000005        1.3785156249999995        1.3607324603828621     
  0.91000000000000003        1.3769531249999998        1.3673715700783504     
   1.0100000000000000        1.3762499999999998        1.3711196345637411     
```

## Rayleigh wave inversion

in preparation

---

## Install
Type `make` in the `src` directory.
