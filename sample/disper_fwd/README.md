# Dispersion curve forward computation


## How to run

`../../bin/disper_fwd disper_fwd.in`

## Inputs

### Main parameter file

* The filename must be specified as the first command line argument.
* The following can be used as a template.

```
#-----------------------------------------------------------------------
# Sample input to disper_fwd of SEIS_FILO package
#-----------------------------------------------------------------------

#---- Required parameters ----------------------------------------------

# Output filename
disper_out     = disper_fwd.out

# Input velocity model
vmod_in        = ../vmod/vmod.in

# Output frequency or period range 
freq_or_period = freq    # Output observables (Freq or Period)
xmin           = 0.01    # Lower bound [Hz] / [s]
xmax           = 1.0     # Upper bound [Hz] / [s]
dx             = 0.05    # Interval    [Hz] / [s]

# Phase velocity range for root search
cmin           = 1.5     # Lower bound [km/s]
cmax           = 4.6     # Upper bound [km/s]
dc             = 0.02    # Interval    [km/s]

# Phase type
disper_phase   = R       # Currently only 'R' (Rayleigh wave) is supported
n_mode         = 0       # Mode number 
                         # (0: fundamental mode, 1: 1st higher mode, 
                         #  2: 2nd higher mode, ...)

#---- Optional parameters ----------------------------------------------

# Additive noise 
noise_added    = 0.0     # Standard deviation of noise (default=0.0)

# Random number seeds (used for noise generation)
i_seed1        = 1001
i_seed2        = 22002
i_seed3        = 333003
i_seed4        = 4444004
```

### Velocity model file

* The filename must be specified in the parameter file.
* The following can be used as a template.
* For the ocean layer, set Vs to a negative value. This is allowed only for topmost layer.
* Thickness of the bottom layer is ignored.

```
#-----------------------------------------------------------------------
# Velocity model file:
# Sample input to disper_fwd & recv_func_fwd of SEIS_FILO package
#-----------------------------------------------------------------------
# Number of layers (including ocean & bottom layers)
         4
# Vp(km/s) Vs(km/s) Density(g/cm^3) Thickness(km)
       1.5     -3.0             1.0           2.0
       5.0      3.0             2.5           5.0
       7.0      4.0             3.0          10.0
       8.0      4.6             3.3           0.0
```

## Output

### Dispersion curve file

* Each row contains frquency (or period), phase velocity, group velocity, and ellipticity.
* The following is a sample output
```
    0.0100    4.1693    4.1138    0.7792
    0.0600    3.9255    3.6000    0.8886
    0.1100    3.5184    2.6941    0.7256
    0.1600    2.9291    1.6103    0.5378
    0.2100    2.2273    1.1432    0.4644
    0.2600    1.8949    1.1989    0.4441
    0.3100    1.7445    1.2681    0.4315
    0.3600    1.6643    1.3201    0.4241
    0.4100    1.6165    1.3573    0.4197
    0.4600    1.5857    1.3843    0.4169
    0.5100    1.5648    1.4041    0.4151
    0.5600    1.5497    1.4194    0.4138
    0.6100    1.5387    1.4310    0.4129
    0.6600    1.5304    1.4402    0.4122
    0.7100    1.5240    1.4476    0.4117
    0.7600    1.5189    1.4537    0.4113
    0.8100    1.5149    1.4587    0.4110
    0.8600    1.5117    1.4627    0.4107
    0.9100    1.5090    1.4662    0.4105
    0.9600    1.5068    1.4692    0.4103
    1.0100    1.5049    1.4717    0.4102
```
